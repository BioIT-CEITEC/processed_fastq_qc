
def biobloom_input(wildcards):
    preprocessed = "processed_fastq"
    input = {}
    if not config["is_paired"]:
        input['r1'] = os.path.join(preprocessed,"{sample}.fastq.gz")
    else:
        input['r1'] = os.path.join(preprocessed,"{sample}_R1.fastq.gz")
        input['r2'] = os.path.join(preprocessed,"{sample}_R2.fastq.gz")
    return  input

rule biobloom:
    input:  unpack(biobloom_input)
    output: table = "qc_reports/{sample}/biobloom/{sample}.biobloom_summary.tsv",
    log:    "logs/{sample}/biobloom.log",
    threads: 8
    resources: mem=30
    params: prefix = "qc_reports/{sample}/biobloom/{sample}.biobloom",
            ref_list = config["biobloom_ref"],
            ref_dir = GLOBAL_REF_PATH,
            paired = paired,
    conda: "../wrappers/biobloom/env.yaml"
    script: "../wrappers/biobloom/script.py"


rule merge_detected_species:
    input:  table = expand("qc_reports/{sample}/species_detector/{sample}_{read_pair_tag}.species_stats.tsv", sample = sample_tab.sample_name, read_pair_tag = read_pair_tags),
    output: table = "qc_reports/species_detector_summary_mqc.tsv",
    log:    "logs/merge_detected_species.log",
    params: nreads = config["max_reads_for_sp_detector"],
            evalue = config["evalue_for_sp_detector"],
            tmpd = GLOBAL_TMPD_PATH,
            top_sp = 10,
    conda:  "../wrappers/merge_detected_species/env.yaml"
    script: "../wrappers/merge_detected_species/script.R"


def species_detector_input(wcs):
    if wcs.read_pair_tag == "SE":
      return f'processed_fastq/{wcs.sample}.fastq.gz'
    else:
      return f'processed_fastq/{wcs.sample}_{wcs.read_pair_tag}.fastq.gz'

rule species_detector:
    input:  fastq = species_detector_input,
            taxdbd= GLOBAL_REF_PATH+"/general/nucleotide_DB/taxdb.btd",
            taxdbi= GLOBAL_REF_PATH+"/general/nucleotide_DB/taxdb.bti",
    output: table = "qc_reports/{sample}/species_detector/{sample}_{read_pair_tag}.species_stats.tsv",
    log:    "logs/{sample}/species_detector_{read_pair_tag}.log",
    threads: 4
    resources: mem=20
    params: fasta = "qc_reports/{sample}/species_detector/{sample}_{read_pair_tag}.fasta",
            blast = "qc_reports/{sample}/species_detector/{sample}_{read_pair_tag}.blast",
            ntdb = GLOBAL_REF_PATH+"/general/nucleotide_DB/nt_db",
            nreads = config["max_reads_for_sp_detector"],
            evalue = config["evalue_for_sp_detector"],
            tmpd = GLOBAL_TMPD_PATH,
    conda: "../wrappers/species_detector/env.yaml"
    script: "../wrappers/species_detector/script.R"
