ASSEMBLER = config["ASSEMBLER"]
THREADS = config["THREADS"]

db_name = config["diamondDB"]["name"]
db_url = config["diamondDB"]["url"]
zipped_file = db_url.split('/')[-1]
unzipped_file = zipped_file[:-3]

rule download_uniref:
    """
    Download
    """
    output:
        refdb = "reference/" + unzipped_file
    run:
        shell("wget " + db_url + " -P reference")
        shell("gunzip reference/" + zipped_file)

rule diamond_index:
    input:
        refdb = "reference/" + unzipped_file
    output:
        "reference/" + unzipped_file + ".dmnd"
    threads:
        THREADS
    run:
        shell("diamond makedb -p {threads} --in {input.refdb} -d {input.refdb}")

rule diamond_blastp:
    input:
        refdb = "reference/" + unzipped_file + ".dmnd",
        fasta = expand("reference/assembled/{assembler}_out/genes.faa", assembler = ASSEMBLER)
    output:
        expand("reference/assembled/{assembler}_out/genes.{ref}.out", assembler = ASSEMBLER, ref = db_name)
    run:
        shell("diamond blastp -d {input.refdb} -q {input.fasta} -o {output} -p {threads} -k 1 -f 6 qseqid stitle pident length \
        mismatch gapopen qstart qend sstart send evalue bitscore")

rule annotate_faa:
    input:
        expand("reference/assembled/{assembler}_out/genes.faa", assembler = ASSEMBLER),
        expand("reference/assembled/{assembler}_out/genes.{ref}.out", assembler = ASSEMBLER, ref = db_name)
    output:
        expand("reference/assembled/{assembler}_out/genes_annotated.faa", assembler = ASSEMBLER)
    script:
        "../scripts/annotate_fasta.py"

rule annotate_fna:
    input:
        expand("reference/assembled/{assembler}_out/genes.fna", assembler = ASSEMBLER),
        expand("reference/assembled/{assembler}_out/genes.{ref}.out", assembler = ASSEMBLER, ref = db_name)
    output:
        expand("reference/assembled/{assembler}_out/genes_annotated.fna", assembler = ASSEMBLER)
    script:
        "../scripts/annotate_fasta.py"

rule annotate_gff:
    input:
        expand("reference/assembled/{assembler}_out/genes.gff", assembler = ASSEMBLER),
        expand("reference/assembled/{assembler}_out/genes.{ref}.out", assembler = ASSEMBLER, ref = db_name)
    output:
        expand("reference/assembled/{assembler}_out/genes_annotated.gff", assembler = ASSEMBLER)
    script:
        "../scripts/annotate_gff.py"
