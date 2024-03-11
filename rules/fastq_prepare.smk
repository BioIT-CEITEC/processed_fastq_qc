rule fastq_prepare_SE:
    input:  in_filename = expand("raw_fastq/{sample}_R1.fastq.gz",sample = sample_tab.sample_name)
    output: fastq = temp("umi_fastq/{sample}.fastq.gz")
    log:    "logs/{sample}/fastq_prepare_SE.log"
    params: umi = config["UMI"],
    threads:  1
    conda:  "../wrappers/fastq_prepare_SE/env.yaml"
    script: "../wrappers/fastq_prepare_SE/script.py"

rule fastq_prepare_PE:
    input:  in_filename = expand("raw_fastq/{sample}_R1.fastq.gz",sample = sample_tab.sample_name)
    output: R1 = temp("umi_fastq/{sample}_R1.fastq.gz"),
            R2 = temp("umi_fastq/{sample}_R2.fastq.gz")
    log:    "logs/{sample}/fastq_prepare_PE.log"
    params: umi = config["UMI"],
    threads:  1
    conda:  "../wrappers/fastq_prepare_PE/env.yaml"
    script: "../wrappers/fastq_prepare_PE/script.py"