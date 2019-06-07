#!/usr/bin/env python

"""
align_genome.snakefile

Snakemake workflow for downloading genome from NCBI and
RNA-Seq files from SRA
"""

__authors__ = "David Levy-Booth, Parker Lloyd"
__copyright__ = "Copyright 2019, David Levy-Booth"
__email__ = "dlevyboo@mail.ubc.ca"
__license__ = "GPL3"

import os
from pathlib import Path
THREADS = config["threads"]
TRIMMER = config["TRIMMER"]
ALIGNER = config["ALIGNER"]
METHOD = config["METHOD"]


# rule all:
#     """
#     Collect the main outputs of the workflow.
#     """
#     input:
#         # expand("genome/{genome_id}/{genome_id}_genomic_prokka.gbff",
#         # genome_id=config["genomes"][config["genome_id"]]["url"].split("/")[-1]),
#         expand("genome/{genome_id}/{genome_id}_genomic.fna",
#         genome_id=config["genomes"][config["genome_id"]]["url"].split("/")[-1]),
#
#         expand("transcriptome/reads/{sra_id}_1.untrimmed.fastq", sra_id = config["sample_ids"]),
#         expand("transcriptome/reads/{sra_id}_2.untrimmed.fastq", sra_id = config["sample_ids"])
#         # expand(count_out,
#         # method=METHOD, aligner=ALIGNER, trimmer=TRIMMER)


rule download_genome:
    params:
        genome_url = config["genomes"][config["genome_id"]]["url"]
    output:
        gff = "genome/{genome_id}/{genome_id}_genomic.gff",
        gbff = "genome/{genome_id}/{genome_id}_genomic.gbff",
        cd_fna = "genome/{genome_id}/{genome_id}_cds_from_genomic.fna",
        prep_gff = "genome/{genome_id}/{genome_id}_genomic_salmon.gff",
        prep_cd_fna = "genome/{genome_id}/{genome_id}_cds_from_genomic_salmon.fna",
        fna = "genome/{genome_id}/{genome_id}_genomic.fna"
    log:
        "log/genome/{genome_id}.log"
    benchmark:
        "benchmarks/{genome_id}.download.benchmark.txt"
    run:
        shell("rsync --copy-links --recursive --times --verbose "
        "rsync://{params.genome_url} genome --log-file={log}")
        shell("gunzip -r genome")
        shell("sed 's/locus_tag/gene_id/g' {output.gff} > {output.prep_gff}")
        shell("sed 's/.*\[locus_tag=\([^]]*\)\].*/>\1/g' {output.cd_fna} > {output.prep_cd_fna}")


rule get_SRA_by_accession:
    """
    Retrieve a single-read FASTQ file from SRA (Sequence Read Archive) by run accession number.
    max_reads: Maximal number of reads to download for each sample.
    """
    output:
        "transcriptome/reads/{sra_id}_1.fastq",
        "transcriptome/reads/{sra_id}_2.fastq"
    params:
        max_reads = config["max_reads"]
    threads: THREADS
    benchmark:
        "benchmarks/{sra_id}.download.benchmark.txt"
    shell:
        """
        fasterq-dump {wildcards.sra_id} -O transcriptome/reads --split-files -p
        # This clears a cache where SRA Tools reserves a lot of space
        cache-mgr --clear >/dev/null 2>&1
        """

rule rename_reads:
    input:
        "transcriptome/reads/{sra_id}_1.fastq",
        "transcriptome/reads/{sra_id}_2.fastq"
    output:
        "transcriptome/reads/{sra_id}_1.untrimmed.fastq",
        "transcriptome/reads/{sra_id}_2.untrimmed.fastq"
    run:
        shell("rename 's/.fastq/.untrimmed.fastq/' {input}")
