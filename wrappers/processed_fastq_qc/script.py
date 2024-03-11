######################################
# wrapper for rule: processed_fastq_qc
######################################
import subprocess
from os.path import dirname
from os.path import basename
from snakemake.shell import shell
import re
shell.executable("/bin/bash")
log_filename = str(snakemake.log)

f = open(log_filename, 'wt')
f.write("\n##\n## RULE: processed_fastq_qc \n##\n")
f.close()

version = str(subprocess.Popen("conda list ", shell=True, stdout=subprocess.PIPE).communicate()[0], 'utf-8')
f = open(log_filename, 'at')
f.write("## CONDA: "+version+"\n")
f.close()

command = "mkdir -p " + dirname(snakemake.output.html) + " >> "+log_filename+" 2>&1"
f = open(log_filename, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)

if snakemake.params.paired:
    if re.search("_R1.fastq",str(snakemake.input.processed)):
        tags = "_R1"
    else:
        tags = "_R2"
else:
    tags = ""

tag_input_fastq = dirname(snakemake.input.processed) + "/" + snakemake.wildcards.sample + tags + ".fastq.gz"

f = open(log_filename, 'at')
f.write("## tag_input_fastq: "+tag_input_fastq+"\n")
f.close()

command = "fastqc -o " + dirname(snakemake.output.html) + " "+snakemake.params.extra+" --threads "+str(snakemake.threads)+" "+tag_input_fastq+" >> "+log_filename+" 2>&1 "
f = open(log_filename, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)

command = "mv " + dirname(snakemake.output.html) + "/" + basename(tag_input_fastq).replace(".fastq.gz","_fastqc.html") + " " + snakemake.output.html
f = open(log_filename, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)

if tags == "":
    zip_out_basename = "SE_trim_fastqc.zip"
else:
    zip_out_basename = tags.replace("_","") + "_trim_fastqc.zip"

command = "mv " + dirname(snakemake.output.html) + "/" + basename(tag_input_fastq).replace(".fastq.gz","_fastqc.zip") + " " + dirname(snakemake.output.html) + "/" + zip_out_basename
f = open(log_filename, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)

# command = ' sed -r -i "s:<h2>[ ]*Summary[ ]*<\/h2><ul>:&<li><b>Return to <a href=\'../'+snakemake.params.lib_name+'.final_report.html\'>start page<\/a><\/b><\/li>:" '+snakemake.output.html
# f = open(log_filename, 'at')
# f.write("## COMMAND: "+command+"\n")
# f.close()
# shell(command)
