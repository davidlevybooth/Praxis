#!/usr/bin/env python

__authors__ = "David Levy-Booth, Parker Lloyd"
__copyright__ = "Copyright 2019, David Levy-Booth"
__email__ = "dlevyboo@mail.ubc.ca"
__license__ = "GPL3"

import os
import itertools
import pandas as pd
from pathlib import Path

TRIMMER = config["TRIMMER"]
ALIGNER = config["ALIGNER"]
METHOD = config["METHOD"]
ASSEMBLER = config["ASSEMBLER"]

genome_url = config["genome"]["ncbi_url"]
samples = pd.read_csv("samples.tsv", sep="\t")
contrasts = list(itertools.combinations(set(samples["Condition"]), 2))
contrasts = sorted(['_'.join(map(str,sorted(pair))) for pair in contrasts])


# if "salmon" in METHOD and len(METHOD) == 1:
#     count_out = "results/tables/salmon/{trimmer}/counts.tsv"
#     DE_out = count_out.replace("counts", "{contrasts}")
# if "salmon" in METHOD and len(METHOD) > 1:
#     METHOD.remove("salmon")
#     count_out = ["results/tables/salmon/{trimmer}/counts.tsv",
#     "results/tables/{method}/{trimmer}/{aligner}/counts.tsv"]
#     DE_out = [table.replace("counts", "{contrasts}") for table in count_out]


if METHOD == "salmon":
    count_out = "results/tables/salmon/{trimmer}/counts.tsv"
    DE_out = count_out.replace("counts", "{contrasts}")
else:
    count_out = "results/tables/{method}/{trimmer}/{aligner}/counts.tsv"
    DE_out = count_out.replace("counts", "{contrasts}")

rule all:
    """
    Collect the main outputs of the workflow.
    """
    input:
        expand("reference/assembled/{assembler}_out/genes_annotated.faa", assembler = ASSEMBLER),
        expand("reference/assembled/{assembler}_out/genes_annotated.fna", assembler = ASSEMBLER),
        expand("reference/assembled/{assembler}_out/genes_annotated.gff", assembler = ASSEMBLER),
        # expand("transcriptome/qc/fastqc/{trimmer}/{sra_id}_{num}.html", trimmer=TRIMMER, sra_id=config["sample_ids"], num=["1","2"]),
        # expand(count_out, method=METHOD, aligner=ALIGNER, trimmer=TRIMMER),
        # expand(DE_out, method=METHOD, aligner=ALIGNER, trimmer=TRIMMER, contrasts = contrasts)


include: "rules/download.smk"
include: "rules/index_align.smk"
include: "rules/quality_control.smk"
include: "rules/quantify.smk"
include: "rules/deseq2.smk"
include: "rules/annotate.smk"
if not genome_url:
    if ASSEMBLER:
        include: "rules/assemble.smk"
    else:
        raise ValueError("An assembler must be specified if no genome URL is provided.")
