import os
import itertools
import pandas as pd
from pathlib import Path

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
