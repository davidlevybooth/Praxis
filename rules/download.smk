#!/usr/bin/env python

"""
Snakemake workflow for downloading genome from NCBI and
RNA-Seq files from SRA
"""

__authors__ = "David Levy-Booth, Parker Lloyd"
__copyright__ = "Copyright 2019, David Levy-Booth"
__email__ = "dlevyboo@mail.ubc.ca"
__license__ = "GPL3"

import os
from pathlib import Path
THREADS = config["THREADS"]
TRIMMER = config["TRIMMER"]
ALIGNER = config["ALIGNER"]
METHOD = config["METHOD"]

genome_url = config["genomes"][config["genome_id"]]["url"]

if genome_url:
    rule download_genome:
        output:
            gff = "reference/genome/{genome_id}/{genome_id}_genomic.gff",
            gbff = "reference/genome/{genome_id}/{genome_id}_genomic.gbff",
            cd_fna = "reference/genome/{genome_id}/{genome_id}_cds_from_genomic.fna",
            prep_gbff = "reference/genome/{genome_id}/{genome_id}_genomic_prokka.gbff",
            prep_cd_fna = "reference/genome/{genome_id}/{genome_id}_cds_from_genomic_salmon.fna",
            fna = "reference/genome/{genome_id}/{genome_id}_genomic.fna"
        log:
            "log/genome/{genome_id}.log"
        benchmark:
            "benchmarks/{genome_id}.download.benchmark.txt"
        run:
            try:
                shell("rsync --copy-links --recursive --times --verbose "
                "rsync://{genome_url} reference/genome --log-file={log}")
                shell("gunzip -r reference/genome")
                shell("sed 's/product/protein/g' {output.gbff} > {output.prep_gbff}")
                shell("sed -i 's/locus_tag/product/g' {output.prep_gbff}")
                shell("sed -i 's/locus_tag/gene_id/g' {output.gff}")
                shell("sed 's/.*\[locus_tag=\([^]]*\)\].*/>\\1/g' {output.cd_fna} > {output.prep_cd_fna}")
                touch("{output.fna}")
            except Exception as e:
                # print(e.returncode)
                if e.returncode == 23:
                    print("[ERROR]: Genome URL is possibly invalid. Check URL in config file.")
                    raise

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
        fasterq-dump {wildcards.sra_id} -O transcriptome/reads --split-files
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
