ASSEMBLER = config["ASSEMBLER"]

rule download:
    output:
        "uniref90.fasta"
    run:
        shell("wget ftp://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref90/uniref90.fasta.gz")
        shell("gunzip uniref90.fasta.gz")


rule diamond_index:
    input:
        refdb = "uniref90.fasta"
    output:
        "uniref90.fasta.dmnd"
    run:
        shell("diamond makedb -p {threads} --in {input.refdb} -d {input.refdb}")

rule diamond_annotate:
    input:
        refdb = "uniref90.fasta.dmnd",
        fasta = expand("transcriptome/{assembler}_out/genes.faa", assembler = ASSEMBLER)
    output:
        expand("transcriptome/{assembler}_out/genes.uniref90.out", assembler = ASSEMBLER)
    run:
        shell("diamond blastp -d {input.refdb} -q {input.fasta} -o {output} -p {threads} -k 1 -f 6 qseqid stitle pident length \
        mismatch gapopen qstart qend sstart send evalue bitscore")
