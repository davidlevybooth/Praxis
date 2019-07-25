import os
from pathlib import Path
THREADS = config["THREADS"]
TRIMMER = config["TRIMMER"]
ALIGNER = config["ALIGNER"]
METHOD = config["METHOD"]


rule megahit_assemble:
    input:
        expand("transcriptome/reads/{trimmer}/{sra_id}_1.fastq", sra_id=config["sample_ids"], trimmer=TRIMMER),
        expand("transcriptome/reads/{trimmer}/{sra_id}_2.fastq", sra_id=config["sample_ids"], trimmer=TRIMMER)
    params:
        left = lambda wildcards, input: ",".join(input[0]),
        right = lambda wildcards, input: ",".join(input[1]),
        out = "reference/assembled/megahit_out"
    output:
        ref = "reference/assembled/megahit_out/final.contigs.fa"
    threads:
        THREADS
    run:
        shell("megahit -1 {params.left} -2 {params.right} --k-list 25 --no-mercy --bubble-level 0 --prune-level 3 -o {params.out} -m 300000000000 -t {threads}")

rule trinity_reformat_headers:
    input:
        "transcriptome/reads/{sra_id}_{n}.{trimmer}.fastq"
    output:
        "transcriptome/reads/{sra_id}_{n}.{trimmer}.newheaders.fastq"
    run:
        shell("perl -ne 's/SR\S+ (\S+) .+/$1\/{wildcards.n}/; print' {input} > {output}")

rule trinity_assemble:
    input:
        expand("transcriptome/reads/{trimmer}/{sra_id}_1.newheaders.fastq", sra_id=config["sample_ids"], trimmer=config["TRIMMER"]),
        expand("transcriptome/reads/{trimmer}/{sra_id}_2.newheaders.fastq", sra_id=config["sample_ids"], trimmer=config["TRIMMER"])
    params:
        left = lambda wildcards, input: ",".join(input[0]),
        right = lambda wildcards, input: ",".join(input[1]),
        out = "reference/assembled/trinity_out",
    threads:
        THREADS
    output:
        ref = "{params.out}/Trinity.fasta"
    run:
        shell("Trinity --left {params.left} --right {params.right} -o {params.out} --max_memory 30G --CPU {threads}")

def select_reference():
    if ASSEMBLER=="megahit":
        return "reference/assembled/megahit_out/final.contigs.fa"
    elif ASSEMBLER=="trinity":
        return "reference/assembled/trinity_out/Trinity.fasta"

rule prodigal:
    input:
        ref = select_reference()
    output:
        gff = expand("reference/assembled/{assembler}_out/genes.gff", assembler = ASSEMBLER), #mapping file
        faa = expand("reference/assembled/{assembler}_out/genes.faa", assembler = ASSEMBLER), #protein file * this is what we'll need for annotation
        fna = expand("reference/assembled/{assembler}_out/genes.fna", assembler = ASSEMBLER)  #nucleotide file* this is what we'll need for alignment/quantification
    run:
        shell("prodigal -i {input.ref} -f gff -o {output.gff} -a {output.faa} -d {output.fna}")
