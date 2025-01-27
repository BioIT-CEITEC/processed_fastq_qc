rule fastq_prepare_SE:
    input:  in_filename = "raw_fastq/{sample}_R1.fastq.gz"
    output: fastq = temp("umi_fastq/{sample}.fastq.gz")
    log:    "logs/{sample}/fastq_prepare_SE.log"
    params: config = config,
    threads:  1
    script: "../wrappers/fastq_prepare_SE/script.py"

rule fastq_prepare_PE:
    input:  in_filename = "raw_fastq/{sample}_R1.fastq.gz"
    output: R1 = temp("umi_fastq/{sample}_R1.fastq.gz"),
            R2 = temp("umi_fastq/{sample}_R2.fastq.gz")
    log:    "logs/{sample}/fastq_prepare_PE.log"
    params: config = config,
    threads:  1
    script: "../wrappers/fastq_prepare_PE/script.py"