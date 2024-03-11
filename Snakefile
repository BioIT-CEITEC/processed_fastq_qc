import os
import pandas as pd
import json
from snakemake.utils import min_version

min_version("5.18.0")
configfile: "config.json"

GLOBAL_REF_PATH = config["globalResources"]
GLOBAL_TMPD_PATH = config["globalTmpdPath"]

os.makedirs(GLOBAL_TMPD_PATH, exist_ok=True)

##### BioRoot utilities #####
module BR:
    snakefile: gitlab("bioroots/bioroots_utilities", path="bioroots_utilities.smk",branch="master")
    config: config

use rule * from BR as other_*

##### Config processing #####

sample_tab = BR.load_sample()

read_pair_tags = BR.set_read_pair_qc_tags() # ["SE"] / ["R1", "R2"]
pair_tag = BR.set_read_pair_tags() # [""] / ["_R1", "_R2"]
paired = BR.set_paired_tags() # "SE" / "PE"

wildcard_constraints:
    sample = "|".join(sample_tab.sample_name)

##### Target rules #####
rule all:
    input: "qc_reports/processed_fastq_multiqc.html"

##### Modules #####

include: "rules/trimming.smk"
include: "rules/fastq_prepare.smk"
