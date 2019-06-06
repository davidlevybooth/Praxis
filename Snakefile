#!/usr/bin/env python

"""
master.snakefile
"""

__authors__ = "David Levy-Booth, Parker Lloyd"
__copyright__ = "Copyright 2019, David Levy-Booth"
__email__ = "dlevyboo@mail.ubc.ca"
__license__ = "GPL3"

include: 'snakefiles/download.snakefile'

rule all:
    """
    Collect the main outputs of the workflow.
    """
    input:
        # expand("genome/{genome_id}/{genome_id}_genomic_prokka.gbff",
        # genome_id=config["genomes"][config["genome_id"]]["url"].split("/")[-1]),
        expand("genome/{genome_id}/{genome_id}_genomic.fna.gz",
        genome_id=config["genomes"][config["genome_id"]]["url"].split("/")[-1]),

        # expand("transcriptome/reads/{sra_id}_1.untrimmed.fastq", sra_id = config["sample_ids"]),
        # expand("transcriptome/reads/{sra_id}_2.untrimmed.fastq", sra_id = config["sample_ids"])
        # expand(count_out,
        # method=METHOD, aligner=ALIGNER, trimmer=TRIMMER)

include: "download.snakefile"
include: "index_align.snakefile"
include: "quality_control.snakefile"
include: "quantify.snakefile"
#include: "deseq2.snakefile"
