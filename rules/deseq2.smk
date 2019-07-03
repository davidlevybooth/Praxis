import os
import itertools
import pandas as pd
from pathlib import Path
THREADS = config["THREADS"]
TRIMMER = config["TRIMMER"]
ALIGNER = config["ALIGNER"]
METHOD = config["METHOD"]

samples = pd.read_csv("samples.tsv", sep="\t")
contrasts = list(itertools.combinations(set(samples["Condition"]), 2))
contrasts = sorted(['_'.join(map(str,sorted(pair))) for pair in contrasts])

if "salmon" in METHOD and len(METHOD) == 1:
    count_out = "results/tables/salmon.{trimmer}.counts.tsv"
    DE_out = count_out.replace("counts", "{contrasts}")
if "salmon" in METHOD and len(METHOD) > 1:
    METHOD.remove("salmon")
    count_out = ["results/tables/salmon.{trimmer}.counts.tsv",
    "results/tables/{method}.{aligner}.{trimmer}.counts.tsv"]
    DE_out = [table.replace("counts", "{contrasts}") for table in count_out]

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
    #log:
    #    "log/deseq2.log"
    threads: 1
    script:
        "../scripts/deseq.R"
