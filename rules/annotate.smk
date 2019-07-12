ASSEMBLER = config["ASSEMBLER"]

rule download:
    output:
        "reference/uniref90.fasta"
    run:
        shell("wget ftp://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref90/uniref90.fasta.gz")
        shell("gunzip reference/uniref90.fasta.gz")


rule diamond_index:
    input:
        refdb = "uniref90.fasta"
    output:
        "reference/uniref90.fasta.dmnd"
    run:
        shell("diamond makedb -p {threads} --in {input.refdb} -d {input.refdb}")

rule diamond_blastp:
    input:
        refdb = "reference/uniref90.fasta.dmnd",
        fasta = expand("reference/assembled/{assembler}_out/genes.faa", assembler = ASSEMBLER)
    output:
        expand("reference/assembled/{assembler}_out/genes.uniref90.out", assembler = ASSEMBLER)
    run:
        shell("diamond blastp -d {input.refdb} -q {input.fasta} -o {output} -p {threads} -k 1 -f 6 qseqid stitle pident length \
        mismatch gapopen qstart qend sstart send evalue bitscore")

rule annotate_faa:
    input:
        expand("reference/assembled/{assembler}_out/genes.faa", assembler = ASSEMBLER),
        expand("reference/assembled/{assembler}_out/genes.uniref90.out", assembler = ASSEMBLER)
    output:
        expand("reference/assembled/{assembler}_out/genes_annotated.faa", assembler = ASSEMBLER)
    script:
        "../scripts/annotate_fasta.py"

rule annotate_fna:
    input:
        expand("reference/assembled/{assembler}_out/genes.fna", assembler = ASSEMBLER),
        expand("reference/assembled/{assembler}_out/genes.uniref90.out", assembler = ASSEMBLER)
    output:
        expand("reference/assembled/{assembler}_out/genes_annotated.fna", assembler = ASSEMBLER)
    script:
        "../scripts/annotate_fasta.py"

rule annotate_gff:
    input:
        expand("reference/assembled/{assembler}_out/genes.gff", assembler = ASSEMBLER),
        expand("reference/assembled/{assembler}_out/genes.uniref90.out", assembler = ASSEMBLER)
    output:
        expand("reference/assembled/{assembler}_out/genes_annotated.gff", assembler = ASSEMBLER)
    script:
        "../scripts/annotate_gff.py"
