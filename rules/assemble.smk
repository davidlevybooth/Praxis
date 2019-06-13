
rule megahit_assemble:
    input:
        read1 = "transcriptome/reads/{sra_id}_1.{trimmer}.fastq",
        read2 = "transcriptome/reads/{sra_id}_2.{trimmer}.fastq"
    output:


rule trinity_reformat_headers:
    input:
        "{sra_id}_{n}.trimmomatic.fastq"
    output:
        "{sra_id}_{n}.trimmomatic.newheaders.fastq"
    shell:
        """
        perl -ne 's/SR\S+ (\S+) .+/$1\/{wildcards.n}/; print' {input} > {output}
        """

rule trinity_assemble:
    input:
        read1 = "transcriptome/reads/{sra_id}_1.{trimmer}.fastq",
        read2 = "transcriptome/reads/{sra_id}_2.{trimmer}.fastq"
    output:

    shell:
        """

        """
