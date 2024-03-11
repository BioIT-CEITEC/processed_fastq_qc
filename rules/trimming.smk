rule merge_fastq_qc:
    input: html=expand("qc_reports/{sample}/processed_fastqc/{read_pair_tag}_trim_fastqc.html",sample=sample_tab.sample_name,read_pair_tag=read_pair_tags)
    output: html="qc_reports/processed_fastq_multiqc.html"
    log: "logs/merge_fastq_qc.log"
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


def preprocess_fastq_input(wildcards):
    if config["UMI"] == "no_umi":
        return expand("raw_fastq/{{sample}}{read_tags}.fastq.gz",read_tags=pair_tag)
    else:
        return expand("umi_fastq/{{sample}}{read_tags}.fastq.gz",read_tags=pair_tag)


rule preprocess:
    input: fastq=preprocess_fastq_input,
    output: processed=expand("processed_fastq/{{sample}}{read_tags}.fastq.gz",read_tags=pair_tag),
    log: "logs/{sample}/preprocessing.log"
    threads: 10
    resources: mem=10
    params: trim_adapters=config["trim_adapters"],
        trim_adapter_select=config["trim_adapter_select"],
        adapter_seq=config["adapter_seq"],
        adapter_type=config["adapter_type"],
        max_error=config["max_error"],
        min_overlap=config["min_overlap"],
        quality_trim=config["quality_trim"],
        r1u="processed_fastq/trimmed/{sample}_R1.discarded.fastq.gz",
        r2u="processed_fastq/trimmed/{sample}_R2.discarded.fastq.gz",
        cut_left1=config[
            "cut_left1"],# Applied only if trim left is true, trimming from R1 (different for classic:0, quant:10, sense:9)
        cut_right1=config[
            "cut_right1"],# Applied only if trim right is true, trimming from R1; you should allow this if you want to trim the last extra base and TRIM_LE is true as RD_LENGTH is not effective
        cut_left2=config[
            "cut_left2"],# Applied only if trim left is true, trimming from R2 (different for classic:0, quant:?, sense:7)
        cut_right2=config[
            "cut_right2"],# Applied only if trim right is true, trimming from R2; you should allow this if you want to trim the last extra base and TRIM_LE is true as RD_LENGTH is not effective
        quality_base=config["quality_base"],
        min_length=config["min_length"],
        max_length=config["max_length"],
        trim_stats="qc_reports/{sample}/cutadapt/{sample}_preprocessing.log"
    conda: "../wrappers/preprocess/env.yaml"
    script: "../wrappers/preprocess/script.py"

