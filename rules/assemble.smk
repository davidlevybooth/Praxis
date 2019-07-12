import os
from pathlib import Path
THREADS = config["THREADS"]
TRIMMER = config["TRIMMER"]
ALIGNER = config["ALIGNER"]
METHOD = config["METHOD"]


rule megahit_assemble:
    input:
        expand("transcriptome/reads/{sra_id}_1.{trimmer}.fastq", sra_id=config["sample_ids"], trimmer=config["TRIMMER"]),
        expand("transcriptome/reads/{sra_id}_2.{trimmer}.fastq", sra_id=config["sample_ids"], trimmer=config["TRIMMER"])
    params:
        left = lambda input: ",".join(input[0]),
        right = lambda input: ",".join(input[1])
    output:
        directory("reference/assembled/megahit_out")
    shell:
        """
        megahit -1 {params.left} -2 {params.right} --k-list 25 --no-mercy --bubble-level 0 --prune-level 3 -o {output} -m 300000000000
        """

rule trinity_reformat_headers:
    input:
        "transcriptome/reads/{sra_id}_{n}.trimmomatic.fastq"
    output:
        "transcriptome/reads/{sra_id}_{n}.trimmomatic.newheaders.fastq"
    shell:
        """
        perl -ne 's/SR\S+ (\S+) .+/$1\/{wildcards.n}/; print' {input} > {output}
        """

rule trinity_assemble:
    input:
        expand("transcriptome/reads/{sra_id}_1.{trimmer}.newheaders.fastq", sra_id=config["sample_ids"], trimmer=config["TRIMMER"]),
        expand("transcriptome/reads/{sra_id}_2.{trimmer}.newheaders.fastq", sra_id=config["sample_ids"], trimmer=config["TRIMMER"])
    params:
        left = lambda wildcards, input: ",".join(input[0]),
        right = lambda wildcards, input: ",".join(input[1])
    output:
        directory("reference/assembled/trinity_out")
    shell:
        """
        Trinity --left {params.left} --right {params.right} -o {output} --max_memory 30G --CPU 6
        """

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
    shell:
        """
        prodigal -i {input.ref} -f gff -o {output.gff} -a {output.faa} -d {output.fna}
        sed 's/locus_tag/gene_id/g' {output.gff}
        """
