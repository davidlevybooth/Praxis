#!/usr/bin/env python

__authors__ = "David Levy-Booth, Parker Lloyd"
__copyright__ = "Copyright 2019, David Levy-Booth"
__email__ = "dlevyboo@mail.ubc.ca"
__license__ = "GPL3"

import os
import itertools
import pandas as pd
from pathlib import Path

THREADS = config["THREADS"]
TRIMMER = config["TRIMMER"]
ALIGNER = config["ALIGNER"]
METHOD = config["METHOD"]
ASSEMBLER = config["ASSEMBLER"]

samples = pd.read_csv("samples.tsv", sep="\t")
contrasts = list(itertools.combinations(set(samples["Condition"]), 2))
contrasts = sorted(['_'.join(map(str,sorted(pair))) for pair in contrasts])

include: "rules/download.smk"
include: "rules/index_align.smk"
include: "rules/quality_control.smk"
include: "rules/quantify.smk"
include: "rules/deseq2.smk"
if not ASSEMBLER:
    include: "rules/assemble.smk"

if "salmon" in METHOD and len(METHOD) == 1:
    count_out = "results/tables/salmon.{{trimmer}}.counts.tsv"
    DE_out = count_out.replace("counts", "{contrasts}")
if "salmon" in METHOD and len(METHOD) > 1:
    METHOD.remove("salmon")
    count_out = ["results/tables/salmon.{trimmer}.counts.tsv",
    "results/tables/{method}.{aligner}.{trimmer}.counts.tsv"]
    DE_out = [table.replace("counts", "{contrasts}") for table in count_out]

rule all:
    """
    Collect the main outputs of the workflow.
    """
    input:
        #"results/tables/salmon.trimmomatic.counts.tsv",
        # expand(DE_out, method=METHOD, aligner=ALIGNER, trimmer=TRIMMER, contrasts = contrasts)
        # expand(count_out, method=METHOD, aligner=ALIGNER, trimmer=TRIMMER)
        expand("genome/{genome_id}/{genome_id}_genomic_prokka.gbff",
            genome_id=config["genomes"][config["genome_id"]]["url"].split("/")[-1])
