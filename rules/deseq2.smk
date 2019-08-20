import os
import itertools
import pandas as pd
from pathlib import Path

THREADS = config["THREADS"]
TRIMMER = config["TRIMMER"]
ALIGNER = config["ALIGNER"]
METHOD = config["METHOD"]

samples = pd.read_csv('samples.tsv', sep='\t')
contrasts = list(itertools.combinations(set(samples["Condition"]), 2))
contrasts = sorted(['_'.join(map(str,sorted(pair))) for pair in contrasts])

if METHOD == "salmon":
    count_out = "results/tables/salmon/{trimmer}/counts.tsv"
    DE_out = count_out.replace("counts", "{contrasts}")
else:
    count_out = "results/tables/{method}/{trimmer}/{aligner}/counts.tsv"
    DE_out = count_out.replace("counts", "{contrasts}")

rule deseq2:
    input:
        counts = expand(count_out, method=METHOD, aligner=ALIGNER, trimmer=TRIMMER)
    output:
        tables = expand(DE_out, method=METHOD, aligner=ALIGNER, trimmer=TRIMMER, contrasts = contrasts)
    params:
        samples = samples["SRA"],
        data = config["samples"],
        contrasts = contrasts
    conda:
        "/envs/deseq2.yaml"
    script:
        "../scripts/deseq.R"
