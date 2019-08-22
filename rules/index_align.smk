import os
from pathlib import Path
THREADS = config["THREADS"]
TRIMMER = config["TRIMMER"]
ALIGNER = config["ALIGNER"]
METHOD = config["METHOD"]
ASSEMBLER = config["ASSEMBLER"]

genome_url = config["genome"]["ncbi_url"]
genome_file = config["genome"]["ref_file"]

# Select genome/transcriptome reference directories
if genome_file:
    reference = genome_file
    indexBase = "intermediate/genome_provided"
else:
    if genome_url:
        reference = expand("reference/genome/{genome_id}/{genome_id}_genomic.fna",
            genome_id = genome_url.split("/")[-1])
        indexBase = "intermediate/genome_downloaded"
    else:
        if ASSEMBLER=="megahit":
            reference = "reference/assembled/megahit_out/final.contigs.fa"
            indexBase = "intermediate/megahit/transcriptome"
        elif ASSEMBLER=="trinity":
            reference = "reference/assembled/trinity_out/Trinity.fasta"
            indexBase = "intermediate/trinity/transcriptome"

if "bt2" in ALIGNER:
    rule bt2_index:
        """
        Index a genome using Bowtie 2.
        """
        input:
            ref=reference
        output:
            expand("{indexBase}.{n}.bt2", indexBase = indexBase, n = ["1","2","3","4"]),
            expand("{indexBase}.rev.{n}.bt2", indexBase = indexBase, n = ["1","2"])
        params:
            base = indexBase
        log:
            "log/index.log"
        benchmark:
            "benchmarks/index.bt2.index.benchmark.txt"
        shell:
            """
            bowtie2-build {input} {params.base} 2> {log}
            """

    rule bt2_align:
        """
        Align a fastq file to a genome index using Bowtie 2.
        """
        input:
            fastq_1 = "transcriptome/reads/{trimmer}/{sra_id}_1.fastq",
            fastq_2 = "transcriptome/reads/{trimmer}/{sra_id}_2.fastq",
            index = expand("{indexBase}.{n}.bt2", indexBase = indexBase, n = ["1","2","3","4"]),
            index_rev = expand("{indexBase}.rev.{n}.bt2", indexBase = indexBase, n = ["1","2"])
        output:
            sam = "intermediate/{trimmer}/bt2/{sra_id}.sam",
            bam = "intermediate/{trimmer}/bt2/{sra_id}.bam"
        wildcard_constraints:
            sra_id = '((?!\.).)*'
        params:
            base = indexBase
        threads: THREADS
        log:
            "log/{trimmer}/bt2/{sra_id}.align.log"
        benchmark:
            "benchmarks/{trimmer}/bt2/{sra_id}.align.benchmark.txt"
        run:
            shell("bowtie2 -x {params.base} --threads {threads} -1 {input.fastq_1} -2 {input.fastq_2} -S {output.sam} 2> {log}")
            shell("samtools view -bS {output.sam} > {output.bam}")

if "bbmap" in ALIGNER:
    rule bbmap_align:
        """
        Align a fastq file to a genome index using bbmap.
        """
        input:
            fastq_1 = "transcriptome/reads/{trimmer}/{sra_id}_1.fastq",
            fastq_2 = "transcriptome/reads/{trimmer}/{sra_id}_2.fastq",
            ref = reference
        output:
            sam = "intermediate/{trimmer}/bbmap/{sra_id}.sam",
            bam = "intermediate/{trimmer}/bbmap/{sra_id}.bam"
        threads: THREADS
        log:
            "log/{trimmer}/{sra_id}.bbmap.align.log"
        benchmark:
            "benchmarks/{trimmer}/{sra_id}.bbmap.align.benchmark.txt"
        run:
            shell("bbmap.sh -Xmx20g trimreaddescriptions=t threads={threads} in1={input.fastq_1} in2={input.fastq_2} \
            out={output.sam} ref={input.ref}  path=reference/bbmap_index/ 2> {log}")
            shell("samtools view -bS {output.sam} > {output.bam}")

rule sort_bam:
    """
    Sort a bam file.
    """
    input:
        expand("intermediate/{trimmer}/{aligner}/{{sra_id}}.bam", trimmer=TRIMMER, aligner=ALIGNER),
    output:
        expand("intermediate/{trimmer}/{aligner}/{{sra_id}}.sorted.bam", trimmer=TRIMMER, aligner=ALIGNER),
    threads: THREADS
    shell:
        """
        samtools sort {input} > {output}
        """
