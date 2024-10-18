def merge_fastq_qc_input(wcs):
    inputs = {'html': expand("qc_reports/{sample}/processed_fastqc/{read_pair_tag}_trim_fastqc.html",sample = sample_tab.sample_name, read_pair_tag = read_pair_tags)}
    if config['biobloom']:
        inputs['biobloom'] = expand("qc_reports/{sample}/biobloom/{sample}.biobloom_summary.tsv",sample = sample_tab.sample_name)
    if config['species_detector']:
        inputs['sp_det'] = "qc_reports/species_detector_summary_mqc.tsv"
    return inputs

rule merge_fastq_qc:
    input: unpack(merge_fastq_qc_input)
    output: html="qc_reports/processed_fastq_multiqc.html"
    log: "logs/merge_fastq_qc.log"
    params: trim_adapters = config["trim_adapters"]
    conda: "../wrappers/merge_fastq_qc/env.yaml"
    script: "../wrappers/merge_fastq_qc/script.py"


def processed_fastq_qc_input(wildcards):
    preprocessed = "processed_fastq"
    if read_pair_tags == ["SE"]:
        return os.path.join(preprocessed,"{sample}.fastq.gz")
    else:
        return os.path.join(preprocessed,"{sample}_{read_pair_tags}.fastq.gz")


rule processed_fastq_qc:
    input: processed=processed_fastq_qc_input,
    output: html="qc_reports/{sample}/processed_fastqc/{read_pair_tags}_trim_fastqc.html",
    log: "logs/{sample}/processed_fastqc_{read_pair_tags}.log"
    params: extra="--noextract --format fastq --nogroup",
        paired=config["is_paired"]
    threads: 2
    conda: "../wrappers/processed_fastq_qc/env.yaml"
    script: "../wrappers/processed_fastq_qc/script.py"
