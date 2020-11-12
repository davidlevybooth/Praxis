#!/usr/bin/python3


# Module Imports
#----------------------------------------------------------------------------#

include: "rules/common.smk"

import os
from os import listdir
from os.path import isfile, join
import itertools
import pandas as pd
from pathlib import Path
import argparse
import psutil
import copy
from datetime import datetime
now = datetime.now().strftime("%d/%m/%Y %H:%M:%S")

# Authorship
#----------------------------------------------------------------------------#

__authors__ = "David Levy-Booth, Parker Lloyd"
__copyright__ = "Copyright 2019, David Levy-Booth"
__email__ = "dlevyboo@mail.ubc.ca"
__license__ = "GPL3"


# Configuration 
#----------------------------------------------------------------------------#

TRIMMER = config["TRIMMER"]
ALIGNER = config["ALIGNER"]
METHOD = config["METHOD"]
ASSEMBLER = config["ASSEMBLER"]
THREADS = config["THREADS"]

if config["genome"]["ncbi_url"]:
    genome_url = config["genome"]["ncbi_url"]

if config["genome"]["ref_file"]:
    genome_file = config["genome"]["ref_file"]

samples_table = pd.read_csv(config["samples"], sep="\t")
sampleID_list = samples_table["SampleID"]

sample_directory = config["sample_directory"]
output_directory = config["output_directory"]

contrasts = get_contrasts(samples_table["Condition"])

# Obtain the correct count table for the differential expression analysis
if METHOD == "salmon":
    count_out = output_directory + "results/tables/salmon/{trimmer}/counts.tsv"
    DE_out = count_out.replace("counts", "{contrasts}")
else:
    count_out = output_directory + "results/tables/{method}/{trimmer}/{aligner}/counts.tsv"
    DE_out = count_out.replace("counts", "{contrasts}")


# Begin Snakemake 
#----------------------------------------------------------------------------#

rule all:
    """
    Collect the main outputs of the workflow.
    """
    input:
        # expand("transcriptome/qc/fastqc/{trimmer}/{sra_id}_{num}.html", trimmer=TRIMMER, sra_id=config["sample_ids"], num=["1","2"]),
        expand(count_out, method=METHOD, aligner=ALIGNER, trimmer=TRIMMER),
        expand(DE_out, method=METHOD, aligner=ALIGNER, trimmer=TRIMMER, contrasts = contrasts)

include: "rules/download.smk"
include: "rules/index_align.smk"
include: "rules/quality_control.smk"
include: "rules/quantify.smk"
include: "rules/deseq2.smk"
include: "rules/annotate.smk"
if not genome_url or genome_file:
    if ASSEMBLER:
        include: "rules/assemble.smk"
    else:
        raise ValueError("An assembler must be specified if no genome URL or file is provided.")
