#!/usr/bin/env python

__authors__ = "David Levy-Booth, Parker Lloyd"
__copyright__ = "Copyright 2019, David Levy-Booth"
__email__ = "dlevyboo@mail.ubc.ca"
__license__ = "GPL3"

import os
from pathlib import Path

THREADS = config["threads"]
TRIMMER = config["TRIMMER"]
ALIGNER = config["ALIGNER"]
METHOD = config["METHOD"]

if "salmon" in METHOD and len(METHOD) == 1:
    count_out = "results/tables/salmon.{{trimmer}}.counts.tsv"
if "salmon" in METHOD and len(METHOD) > 1:
    METHOD.remove("salmon")
    count_out = ["results/tables/salmon.{trimmer}.counts.tsv",
    "results/tables/{method}.{aligner}.{trimmer}.counts.tsv"]


rule all:
    """
    Collect the main outputs of the workflow.
    """
    input:
        "results/tables/salmon.trimmomatic.counts.tsv"
        #expand(count_out, method=METHOD, aligner=ALIGNER, trimmer=TRIMMER)

include: "rules/download.smk"
include: "rules/index_align.smk"
include: "rules/quality_control.smk"
include: "rules/quantify.smk"
#include: "deseq2.smk"
