ASSEMBLER = config["ASSEMBLER"]
THREADS = config["THREADS"]

db_name = config["diamondDB"]["name"]
db_url = config["diamondDB"]["url"]
zipped_file = db_url.split('/')[-1]
unzipped_file = zipped_file[:-3]

genome_file = config["genome"]["ref_file"]

# Obtain the correct directory depending on the selected assembler
if genome_file:
    ref = genome_file
    dir = "/".join(genome_file.split("/")[:-2])
else:
    if ASSEMBLER=="megahit":
        ref = "final.contigs.fa"
        dir = "reference/assembled/megahit_genes"
    elif ASSEMBLER=="trinity":
        ref = "Trinity.fasta"
        dir = "reference/assembled/trinity_genes"

rule prodigal:
    """
    Obtain the predicted genes/proteins from the assembled transcriptome.
    """
    input:
        assembly = expand("reference/assembled/{assembler}_out", assembler = ASSEMBLER)
    params:
        reference = "{input.assembly}/" + ref
    output:
        gff = expand("{directory}/genes.gff", directory = dir), #mapping file
        faa = expand("{directory}/genes.faa", directory = dir), #protein file
        fna = expand("{directory}/genes.fna", directory = dir)  #nucleotide file
    run:
        shell("prodigal -i {params.reference} -f gff -o {output.gff} -a {output.faa} -d {output.fna}")

rule download_uniref:
    """
    Download uniref database
    """
    output:
        refdb = "reference/" + unzipped_file
    run:
        shell("wget " + db_url + " -P reference")
        shell("gunzip reference/" + zipped_file)

rule diamond_index:
    """
    Index uniref database
    """
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
        fasta = expand("{directory}/genes.faa", directory = dir)
    output:
        expand("{directory}/genes.{ref}.out", directory = dir, ref = db_name)
    run:
        shell("diamond blastp -d {input.refdb} -q {input.fasta} -o {output} -p {threads} -k 1 -f 6 qseqid stitle pident length \
        mismatch gapopen qstart qend sstart send evalue bitscore")

rule annotate_faa:
    """
    Adding annotations to the fasta headers
    """
    input:
        expand("{directory}/genes.faa", directory = dir),
        expand("{directory}/genes.{ref}.out", directory = dir, ref = db_name)
    output:
        expand("{directory}/genes_annotated.faa", directory = dir)
    script:
        "../scripts/annotate_fasta.py"

rule annotate_fna:
    """
    Adding annotations to the fasta headers
    """
    input:
        expand("{directory}/genes.fna", directory = dir),
        expand("{directory}/genes.{ref}.out", directory = dir, ref = db_name)
    output:
        expand("{directory}/genes_annotated.fna", directory = dir)
    script:
        "../scripts/annotate_fasta.py"

rule annotate_gff:
    """
    Adding annotations to the gff headers
    """
    input:
        expand("{directory}/genes.gff", directory = dir),
        expand("{directory}/genes.{ref}.out", directory = dir, ref = db_name)
    output:
        expand("{directory}/genes_annotated.gff", directory = dir)
    script:
        "../scripts/annotate_gff.py"
