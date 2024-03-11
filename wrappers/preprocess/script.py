######################################
# wrapper for rule: preprocess
######################################
import subprocess
import os
from snakemake.shell import shell
shell.executable("/bin/bash")
log_filename = str(snakemake.log)

f = open(log_filename, 'wt')
f.write("\n##\n## RULE: preprocess \n##\n")
f.close()

version = str(subprocess.Popen("conda list ", shell=True, stdout=subprocess.PIPE).communicate()[0], 'utf-8')
f = open(log_filename, 'at')
f.write("## CONDA: "+version+"\n")
f.close()

command = "mkdir -p " + os.path.dirname(snakemake.params.trim_stats)
f = open(log_filename, 'at')
f.write("## COMMAND: " + command + "\n")
f.close()
shell(command)

if len(snakemake.input.fastq) == 2:
    is_paired = True
    fastq_r1 = snakemake.input.fastq[0]
    fastq_r2 = snakemake.input.fastq[1]
else:
    is_paired = False
    fastq_r1 = snakemake.input.fastq[0]
    fastq_r2 = ""

if len(snakemake.output.processed) == 2:
    fastq_c1 = snakemake.output.processed[0]
    fastq_c2 = " -p " + snakemake.output.processed[1]
    fastq_u1 = snakemake.params.r1u
    fastq_u2 = " --too-short-paired-output " + snakemake.params.r2u
else:
    fastq_c1 = snakemake.output.processed[0]
    fastq_c2 = ""
    fastq_u1 = snakemake.params.r1u
    fastq_u2 = ""

command = "mkdir -p " + os.path.dirname(snakemake.params.r1u)
f = open(log_filename, 'at')
f.write("## COMMAND: " + command + "\n")
f.close()
shell(command)

# Set the command part related to cut the reads
cut_flags = ""
if int(snakemake.params.cut_left1) != 0 or int(snakemake.params.cut_right1) != 0:
    cut_flags = " -u "+str(abs(snakemake.params.cut_left1)) if int(snakemake.params.cut_left1) != 0 else ""
    cut_flags += " -u "+str(snakemake.params.cut_right1) if int(snakemake.params.cut_right1) != 0 else ""

if is_paired:
    if int(snakemake.params.cut_left1) != 0 or int(snakemake.params.cut_right1) != 0 or int(snakemake.params.cut_left2) != 0 or int(snakemake.params.cut_right2) != 0:
        cut_flags = " -u " + str(abs(snakemake.params.cut_left1)) if int(snakemake.params.cut_left1) != 0 else ""
        cut_flags += " -u " + str(-abs(snakemake.params.cut_right1)) if int(snakemake.params.cut_right1) != 0 else ""
        cut_flags += " -U " + str(abs(snakemake.params.cut_left2)) if int(snakemake.params.cut_left2) != 0 else ""
        cut_flags += " -U " + str(-abs(snakemake.params.cut_right2)) if int(snakemake.params.cut_right2) != 0 else ""

simpleClipThreshold = 10
# TODO: check for better settings (see: http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/TrimmomaticManual_V0.32.pdf starting at page 5, or http://www.usadellab.org/cms/?page=trimmomatic)

if snakemake.params.trim_adapters:
  if snakemake.params.trim_adapter_select == "illumina":
    #adapter_list = "1-AGATCGGAAGAGCACACGTCTGAACTCCAGTCA,2-AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"
    adapter_list = "1-AGATCGGAAGAGCACACGTCT,2-AGATCGGAAGAGCGTCGTGTA"
  if snakemake.params.trim_adapter_select == "nextera":
    adapter_list = "CTGTCTCTTATACACATCT"
  if snakemake.params.trim_adapter_select == "smallRNA":
    adapter_list = "TGGAATTCTCGGGTGCCAAGG"
  if snakemake.params.trim_adapter_select == "custom":
    adapter_list = str(snakemake.params.adapter_seq)
else:
    adapter_list = ""

print(adapter_list)

if adapter_list != "":
    adapters = adapter_list.split(",")
    with open(os.path.dirname(fastq_c1) + "/adapter1.fa", "w") as adapter_file1:
      with open(os.path.dirname(fastq_c1) + "/adapter2.fa", "w") as adapter_file2:
        for i, adapter in enumerate(adapters):
            if adapter.split("-")[0] != "2":
                    adapter_file1.write(">adapt" + str(i) + "\n")
                    adapter_file1.write(adapter.replace("1-","") + "\n")
            if adapter.split("-")[0] != "1":
                adapter_file2.write(">adapt" + str(i) + "\n")
                adapter_file2.write(adapter.replace("2-","") + "\n")
    adapter_flags = " -" + snakemake.params.adapter_type + " file:" + os.path.dirname(fastq_c1) + "/adapter1.fa "
    if is_paired:
      adapter_flags = adapter_flags + " -" + snakemake.params.adapter_type.upper() + " file:" + os.path.dirname(fastq_c1) + "/adapter2.fa "
else:
    adapter_flags = ""

command = "cutadapt -j " + str(snakemake.threads) + " --quality-base=" + str(snakemake.params.quality_base) + " \
                -q " + str(snakemake.params.quality_trim) + " -m " + str(snakemake.params.min_length)+ " \
                --too-short-output " + fastq_u1 + fastq_u2 + "\
                -M " + str(snakemake.params.max_length) + cut_flags + adapter_flags + " \
                -o " + fastq_c1 + fastq_c2 + " " + fastq_r1 + " " + fastq_r2 + " >> " + str(snakemake.params.trim_stats) + " 2>&1 "

f = open(log_filename, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)
