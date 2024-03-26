#########################################################################
######################## U M I   P I P E L I N E ########################
#########################################################################

########################
### FIXED PARAMETERS ###
########################
reference_fasta = config.get("reference_fasta")
if not reference_fasta:
    raise RuntimeError("No reference FASTA found. Please specify 'reference_fasta' in config file")

input_folder = config.get('input_fastq')
if not input_folder:
    raise RuntimeError("No input FASTQ files found. Please specify 'input_fastq' in config file")

target_bed = config.get('targets_bed')
if not target_bed:
    raise RuntimeError("No target BED file found. Please spcify 'targets_bed' in config file")

sample_name = config.get("sample_name", "umi_sample")

#########################
## Optional parameters ##
#########################
allowed_umi_errors = config.get("umi_errors", 3)
min_reads_per_cluster = config.get("min_reads_per_cluster", 20)
max_reads_per_cluster = config.get("max_reads_per_cluster", 60)
min_overlap = config.get("min_overlap", 0.9)
balance_strands = config.get("balance_strands", True)
mm = config.get("medaka_model", "r941_min_high_g360")
fwd_context = config.get("fwd_context", "GTATCGTGTAGAGACTGCGTAGG")
rev_context = config.get("rev_context", "AGTGATCGAGTCAGTGCGAGTG")
fwd_umi = config.get("fwd_umi", "TTTVVVVTTVVVVTTVVVVTTVVVVTTT")
rev_umi = config.get("rev_umi", "AAABBBBAABBBBAABBBBAABBBBAAA")
min_length = config.get("min_length", 40)
max_length = config.get("max_length", 60)
filter_reads = config.get("filter_reads", False)
min_read_len = config.get("min_read_len", 100)
min_mean_qual = config.get("min_mean_qual", 70)
varscan_params = config.get("varscan_params", '--variants 1 --output-vcf 1 --min-coverage 8 --min-avg-qual 0 --min-var-freq 0.01 --strand-filter 0 --p-value 1 --min-reads2 2')

########################
########################
########################

def read_bed_names(filename):
    names = []
    with open(filename) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            cols = line.split("\t")
            if len(cols) < 4:
                print("Warning: ignoring {}. No name found".format(line))
                continue
            names.append(cols[3])
    return names

target = read_bed_names(target_bed)
print("Targets: {}".format(" ".join(target)), file=sys.stderr)

balance_strands_param = "--balance_strands"
if not balance_strands:
    balance_strands_param = ""

########################
######### RULES ########
########################

# Run pipeline
rule variants:
    input:
        expand("{name}/targets.bed", name=sample_name),
        expand("{name}/07_align_consensus/{target}_final.bam.bai", name=sample_name, target=target),
        expand("{name}/09_stats/03_{target}_vsearch_cluster_stats.tsv", name=sample_name, target=target),
        expand("{name}/09_stats/04_{target}_final_cluster_stats.tsv", name=sample_name, target=target),
        expand("{name}/08_variants/{target}_final.vcf", name=sample_name, target=target),
        expand("{name}/10_plots/{target}/final_stats.csv", name=sample_name, target=target)

# Copy BED file to work directory
rule copy_bed:
    input:
        target_bed
    output:
        "{name}/targets.bed"
    shell:
        "cp {input} {output}"

# Filter reads (optional)
rule filter_reads:
    input:
        FQ = input_folder
    params:
        min_read_len = min_read_len,
        min_mean_qual = min_mean_qual,
        filter_reads = filter_reads
    output:
        FQ = "{name}/read.filt.fastq.gz",
        STATS = "{name}/09_stats/01_reads_stats.txt"
    threads: 1
    shell:
        """
        printf 'Total reads in file pre filtering: ' 2>&1 | tee {output.STATS}
        if [[ {input.FQ} =~ \.gz$ ]]
        then
            zcat {input.FQ} | echo $((`wc -l`/4)) 2>&1 | tee -a {output.STATS}
        else
            cat {input.FQ} | echo $((`wc -l`/4)) 2>&1 | tee -a {output.STATS}
        fi
        
        if [[ {params.filter_reads} == "True" ]]
        then
            filtlong --min_length {params.min_read_len} --min_mean_q {params.min_mean_qual} {input.FQ} | gzip > {output.FQ}
            printf 'Total reads in file post filtering: ' 2>&1 | tee -a {output.STATS}
            zcat {output.FQ} | echo $((`wc -l`/4)) 2>&1 | tee -a {output.STATS}
        else
            cp {input.FQ} {output.FQ}
        fi
        """

# Align all reads to reference genome
rule map_1d:
    input:
        FQ = "{name}/read.filt.fastq.gz",
        REF = reference_fasta
    output:
        BAM = "{name}/01_align/1_d.bam",
        BAI = "{name}/01_align/1_d.bam.bai"
    threads: 30
    shell:
        """
        minimap2 -ax map-ont -k 13 -t {threads} {input.REF} {input.FQ} | samtools sort -@ 5 -o {output.BAM} - && samtools index -@ {threads} {output.BAM}
        rm {input.FQ} #because this is a copy now
        """

# Split reads by amplicons
rule split_reads:
    input:
        "{name}/01_align/1_d.bam"
    output:
        DIR = directory("{name}/02_amplicon_fastq/"),
        STATS = "{name}/09_stats/02_umi_filter_reads_stats.txt"
    params:
        bed = target_bed,
        min_overlap = min_overlap
    shell:
        """
        mkdir -p {output.DIR}
        umi_filter_reads --min_overlap {params.min_overlap} -o {output.DIR} {params.bed} {input} 2>&1 | tee {output.STATS}
        """

# Extract UMI sequences from amplicon reads
rule detect_umi_fasta:
    input:
        "{name}/02_amplicon_fastq/"
    output:
        "{name}/03_umi_fasta/{target}_detected_umis.fasta"
    params:
        errors = allowed_umi_errors,
        fwd_context = fwd_context,
        rev_context = rev_context,
        fwd_umi = fwd_umi,
        rev_umi = rev_umi,
    shell:
        """
        umi_extract --fwd-context {params.fwd_context} --rev-context {params.rev_context} --fwd-umi {params.fwd_umi} --rev-umi {params.rev_umi} --max-error {params.errors} {input}/{wildcards.target}.fastq -o {output} --tsv {output}.tsv
        """

# Cluster reads based on UMIs
rule cluster:
    input: 
        "{name}/03_umi_fasta/{target}_detected_umis.fasta"
    output:
        CENT = "{name}/04_umi_clustering/{target}/clusters_centroid.fasta",
        CONS = "{name}/04_umi_clustering/{target}/clusters_consensus.fasta",
        DIR = directory("{name}/04_umi_clustering/{target}/vsearch_clusters")
    params:
        min_length = min_length,
        max_length = max_length
    threads: 10
    shell:
        """
        mkdir -p {wildcards.name}/04_umi_clustering/{wildcards.target}/vsearch_clusters
        vsearch --clusterout_id --clusters {wildcards.name}/04_umi_clustering/{wildcards.target}/vsearch_clusters/test --centroids {output.CENT} --consout {output.CONS} --minseqlength {params.min_length} --maxseqlength {params.max_length} --qmask none --threads {threads} --cluster_fast {input} --clusterout_sort --gapopen 0E/5I --gapext 0E/2I --mismatch -8 --match 6 --iddef 0 --minwordmatches 0 -id 0.85
        """

# Filter clusters based on min/max reads per cluster
rule reformat_filter_clusters:
    input:
        "{name}/04_umi_clustering/{target}/clusters_consensus.fasta",
        "{name}/04_umi_clustering/{target}/vsearch_clusters"
    params:
        min_reads_per_cluster = min_reads_per_cluster,
        max_reads_per_cluster = max_reads_per_cluster,
        balance_strands_param = balance_strands_param
    output:
        out_dir = directory("{name}/05_cluster_filtering/{target}/clusters_fa/"),
        stats = "{name}/09_stats/03_{target}_vsearch_cluster_stats.tsv",
        out_file = "{name}/05_cluster_filtering/{target}/smolecule_clusters.fa"
    shell:
        """
        umi_parse_clusters --smolecule_out {output.out_file} {params.balance_strands_param} --min_reads_per_clusters {params.min_reads_per_cluster} --max_reads_per_clusters {params.max_reads_per_cluster} --stats_out {output.stats} -o {output.out_dir} {input}
        """

# Polish clusters 
rule polish_clusters:
    input:
        "{name}/05_cluster_filtering/{target}/smolecule_clusters.fa"
    output:
        FOLDER = directory("{name}/06_polish_clusters/{target}_consensus_tmp"),
        BAM = "{name}/06_polish_clusters/{target}_consensus.bam",
        F = "{name}/06_polish_clusters/{target}_consensus.fasta"
    params:
        medaka_model = mm
    threads: 30
    shell:
        """
        rm -rf {output.FOLDER}
        medaka smolecule --threads {threads} --length 50 --depth 2 --model {params.medaka_model} --method spoa {output.FOLDER} {input} 2> {output.BAM}_smolecule.log
        mv {output.FOLDER}/consensus.fasta {output.F}
        mv {output.FOLDER}/subreads_to_spoa.bam {output.BAM}
        mv {output.FOLDER}/subreads_to_spoa.bam.bai {output.BAM}.bai
        """

# Reformat fasta header
rule reformat_consensus_clusters_test:
    input:
        "{name}/06_polish_clusters/{target}_consensus.fasta"
    output:
        "{name}/06_polish_clusters/{target}_final.fasta"
    shell:
        """
        sed '/^>/ s/ .*//' {input} > {output}
        """

# Map consensus reads after polishing
rule map_consensus:
    input:
        FA = "{name}/06_polish_clusters/{target}_final.fasta",
        REF = reference_fasta
    output:
        BAM = "{name}/07_align_consensus/{target}_{type}.bam",
        BAI = "{name}/07_align_consensus/{target}_{type}.bam.bai"
    threads: 3
    shell:
        "minimap2 -ax map-ont -k 13 -t {threads} {input.REF} {input.FA} | samtools sort -@ 5 -o {output.BAM} - && samtools index -@ {threads} {output.BAM}"

# Generate stats for final clusters
rule seqkit_bam_acc_tsv:
    input:
        "{name}/07_align_consensus/{target}_final.bam"
    output:
        "{name}/09_stats/04_{target}_final_cluster_stats.tsv"
    shell:
        """
        echo -e "Read\tCluster_size\tChr\tStart\tEnd\tMapQual\tAcc\tReadLen\tRefLen\tRefAln\tRefCov\tReadAln\tReadCov\tStrand\tMeanQual\tLeftClip\tRightClip\tFlags\tIsSec\tIsSup" > {output} && seqkit bam {input} 2>&1 | sed 's/_/\t/' | tail -n +2 >> {output}
        """

# Call variants
rule call_variants:
    input:
        BAM = "{name}/07_align_consensus/{target}_final.bam",
        BED = "{name}/targets.bed",
        REF = reference_fasta
    output:
        "{name}/08_variants/{target}_final.vcf"
    params:
        varscan_params = varscan_params
    shell:
        "samtools mpileup -q 0 -Q 0 -B -d 10000000 -A -l {input.BED} -f {input.REF} {input.BAM} | varscan mpileup2cns {params.varscan_params} > {output}"

# Generate summary stats and plots
rule stats_plots:
    input:
        VCF = "{name}/08_variants/{target}_final.vcf",
    output:
        DIR = directory("{name}/10_plots/{target}"),
        CSV = "{name}/10_plots/{target}/final_stats.csv",
        PDF = "{name}/10_plots/{target}/cluster_size_distribution.pdf"
    params:
        sample_name = sample_name,
        min_reads_per_cluster = min_reads_per_cluster,
        max_reads_per_cluster = max_reads_per_cluster
    shell:
        """
        mkdir -p {wildcards.name}/10_plots/{target}
        umi_stats -m {params.min_reads_per_cluster} -M {params.max_reads_per_cluster} -n {target} -f {output.DIR} {params.sample_name}
        """
        