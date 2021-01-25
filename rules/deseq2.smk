import os
import itertools
import pandas as pd
from pathlib import Path

if config["genome"]["ncbi_url"]:
    genome = config["genome"]["ncbi_url"]
if config["genome"]["ref_file"]:
    genome = config["genome"]["ref_file"]
else:
    genome = "no reference"

if config["project_name"]:
    project_name = config["project_name"]
else:
    project_name = "Differental Expression Analysis"
    
rule deseq2:
    input:
        counts = expand(count_out, method=METHOD, aligner=ALIGNER, trimmer=TRIMMER)
    output:
        tables = expand(DE_out, method=METHOD, aligner=ALIGNER, trimmer=TRIMMER, contrasts = contrasts)
    params:
        samples = samples_table["SampleID"],
        data = config["samples"],
        contrasts = contrasts
    conda:
        "/envs/deseq2.yaml"
    script:
        "../scripts/deseq.R"


rule report:
    input:
        counts = expand(count_out, method=METHOD, aligner=ALIGNER, trimmer=TRIMMER),
        tables = expand(DE_out, method=METHOD, aligner=ALIGNER, trimmer=TRIMMER, contrasts = contrasts),
    output:
        dereport_html = "DE_Report.html"
    params:
        data = config["samples"],
        contrasts = contrasts,
        TRIMMER = config["TRIMMER"],
        ALIGNER = config["ALIGNER"],
        METHOD = config["METHOD"],
        ASSEMBLER = config["ASSEMBLER"],
        GENOME = genome,
        PROJECT = project_name,
    script:
        "../scripts/DE_report.Rmd"
