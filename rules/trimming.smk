# def preprocess_fastq_input(wildcards):
#     if config["UMI"] == "no_umi":
#         return expand("raw_fastq/{{sample}}{read_tags}.fastq.gz",read_tags=pair_tag)
#     else:
#         return expand("umi_fastq/{{sample}}{read_tags}.fastq.gz",read_tags=pair_tag)

rule trimtolibsize:
    input: fastq=expand("umi_fastq/{{sample}}{read_tags}.fastq.gz",read_tags=pair_tag),
    output: processed=expand("libsize_fastq/{{sample}}{read_tags}.fastq.gz",read_tags=pair_tag),
    log: "logs/{sample}/libsizetrimming.log"
    threads: 10
    resources: mem=10
    params: trim_adapters=false, # No adapter trimming
        quality_trim="0", # No quality trimming
        r1u="libsize_fastq/trimmed/{sample}_R1.discarded.fastq.gz",
        r2u="libsize_fastq/trimmed/{sample}_R2.discarded.fastq.gz",
        cut_left1="0",# Applied only if trim left is true, trimming from R1 (different for classic:0, quant:10, sense:9)
        cut_right1=config[
            "lib_forward_read_length"],# Applied only if trim right is true, trimming from R1; you should allow this if you want to trim the last extra base and TRIM_LE is true as RD_LENGTH is not effective
        cut_left2="0",# Applied only if trim left is true, trimming from R2 (different for classic:0, quant:?, sense:7)
        cut_right2=config[
            "lib_reverse_read_length"],# Applied only if trim right is true, trimming from R2; you should allow this if you want to trim the last extra base and TRIM_LE is true as RD_LENGTH is not effective
        quality_base=config["quality_base"],
        min_length="0",
        max_length=config["max_length"],
        trim_stats="qc_reports/{sample}/cutadapt/{sample}_libsizetrimming.log",
        UMI_write_to=config["UMI_write_to"]

    conda: "../wrappers/preprocess/env.yaml"
    script: "../wrappers/preprocess/script.py"

rule preprocess:
    input: fastq=expand("libsize_fastq/{{sample}}{read_tags}.fastq.gz",read_tags=pair_tag),
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
        trim_stats="qc_reports/{sample}/cutadapt/{sample}_preprocessing.log",
        UMI_write_to=config["UMI_write_to"]

    conda: "../wrappers/preprocess/env.yaml"
    script: "../wrappers/preprocess/script.py"

rule filesender:
    input:  processed_fastq = expand("processed_fastq/{sample}{read_tags}.fastq.gz", sample = sample_tab.sample_name,read_tags=pair_tag),
            html = expand("qc_reports/{sample}/processed_fastqc/{read_pair_tags}_trim_fastqc.html",sample=sample_tab.sample_name,read_pair_tags=read_pair_tags)
    output: raw = "raw_fastq.tar.gz",
            processed = "processed_fastq.tar.gz",
            libsize = "libsize_fastq.tar.gz"
    log:    "logs/filesender.log"
    params: recipient = config["recipient"],
            subject = config["entity_name"],
            message = config["message"],
            filesender_raw = config["filesender_raw"],
            filesender_processed = config["filesender_processed"],
            filesender_libsize = config["filesender_libsize"],
            credentials = GLOBAL_REF_PATH + "/reference_info/filesender/filesender_params.json",
            res_file = "qc_reports/",
            res_processed = "processed_fastq/",
            res_libsize = "libsize_fastq/",
            res_raw = "raw_fastq/",
    conda:  "../wrappers/filesender/env.yaml"
    script: "../wrappers/filesender/script.py"
