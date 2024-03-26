"""
UMI amplicion tools setup script

Copyright (c) 2020 by Oxford Nanopore Technologies Ltd.
"""
import os
from setuptools import setup

__version__ = '1.0.0'

setup(
    name='umi-amplicon-tools',
    version=__version__,
    author='Nanopore Technologies Ltd.',
    description='Toolset to work with ONT amplicon sequencing using UMIs',
    zip_safe=False,
    install_requires=[
        'seaborn==0.11.2',
        'tqdm==4.66.2',
        'pysam==0.22.0',
        'pandas==1.4.4',
        'edlib',
        'biopython==1.83',
        'tensorflow==2.7.0',
        'scipy==1.9.1',
        'numpy==1.19.5',
        'matplotlib==3.6'
    ],
    packages=['umi_amplicon_tools'],
    package_data={
        'umi_amplicon_tools': ['data/*'],
    },
    entry_points={
        "console_scripts": [
            'umi_filter_reads = umi_amplicon_tools.filter_reads:main',
            'umi_extract = umi_amplicon_tools.extract_umis:main',
            'umi_parse_clusters = umi_amplicon_tools.parse_clusters:main',
            'umi_bam_to_phred = umi_amplicon_tools.bam_to_phred:main',
            'umi_stats = umi_amplicon_tools.umi_stats:main'
        ]
    },
)

