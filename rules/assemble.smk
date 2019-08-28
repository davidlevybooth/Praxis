import os
from pathlib import Path
THREADS = config["THREADS"]
TRIMMER = config["TRIMMER"]
ALIGNER = config["ALIGNER"]
METHOD = config["METHOD"]

rule megahit_assemble:
    """
    Assemble a reference transcriptome with megahit.
    """
    input:
        left = expand("transcriptome/reads/{trimmer}/{sra_id}_1.fastq", sra_id=config["sample_ids"], trimmer=TRIMMER),
        right = expand("transcriptome/reads/{trimmer}/{sra_id}_2.fastq", sra_id=config["sample_ids"], trimmer=TRIMMER)
    params:
        left = lambda wildcards, input: ",".join(input.left),
        right = lambda wildcards, input: ",".join(input.right),
        memory = config["MAX_MEMORY"]
    output:
        directory("reference/assembled/megahit_out")
    threads:
        THREADS
    run:
        shell("megahit -1 {params.left} -2 {params.right} --k-list 25 --no-mercy --bubble-level 0 --prune-level 3 -m {params.memory} -t {threads} -o {output}")

rule trinity_reformat_headers:
    """
    Remove fasta file header so it can be accepted by trinity.
    """
    input:
        "transcriptome/reads/{sra_id}_{n}.{trimmer}.fastq"
    output:
        "transcriptome/reads/{sra_id}_{n}.{trimmer}.newheaders.fastq"
    run:
        shell("perl -ne 's/SR\S+ (\S+) .+/$1\/{wildcards.n}/; print' {input} > {output}")

rule trinity_assemble:
    """
    Assemble a reference transcriptome with trinity.
    """
    input:
        left = expand("transcriptome/reads/{trimmer}/{sra_id}_1.newheaders.fastq", sra_id=config["sample_ids"], trimmer=config["TRIMMER"]),
        right = expand("transcriptome/reads/{trimmer}/{sra_id}_2.newheaders.fastq", sra_id=config["sample_ids"], trimmer=config["TRIMMER"])
    params:
        left = lambda wildcards, input: ",".join(input.left),
        right = lambda wildcards, input: ",".join(input.right),
        memory = config["MAX_MEMORY"] / 1000000000 # roughly convert to gigabytes
    threads:
        THREADS
    output:
        directory("reference/assembled/trinity_out")
    run:
        shell("Trinity --left {params.left} --right {params.right} -o {output} --max_memory {params.memory} --CPU {threads}")
