######################################
# wrapper for rule: filesender
######################################
import subprocess
import json
from snakemake.shell import shell
shell.executable("/bin/bash")
log_filename = str(snakemake.log)

f = open(log_filename, 'wt')
f.write("\n##\n## RULE: filesender \n##\n")
f.close()

version = str(subprocess.Popen("conda list ", shell=True, stdout=subprocess.PIPE).communicate()[0], 'utf-8')
f = open(log_filename, 'at')
f.write("## CONDA: "+version+"\n")
f.close()

ff = open(str(snakemake.params.credentials))
filesender_credentials = json.load(ff)
ff.close()

username = filesender_credentials["username"]
apikey = filesender_credentials["apikey"]

if snakemake.params.filesender_processed:
    command = "tar -czvf " + snakemake.output.processed + " " + snakemake.params.res_processed + "* >> " + log_filename + " 2>&1"
    f = open(log_filename, 'at')
    f.write("## COMMAND: " + command + "\n")
    f.close()
    shell(command)
    filetosend = " " + snakemake.output.processed
else:
    command = "touch " + snakemake.output.processed + " >> " + log_filename + " 2>&1"
    f = open(log_filename, 'at')
    f.write("## COMMAND: " + command + "\n")
    f.close()
    shell(command)
    filetosend = ""

if snakemake.params.filesender_raw:
    command = "tar -czvf " + snakemake.output.raw + " " + snakemake.params.res_raw + "* >> " + log_filename + " 2>&1"
    f = open(log_filename, 'at')
    f.write("## COMMAND: " + command + "\n")
    f.close()
    shell(command)
    filetosend = " " + snakemake.output.raw
else:
    command = "touch " + snakemake.output.raw + " >> " + log_filename + " 2>&1"
    f = open(log_filename, 'at')
    f.write("## COMMAND: " + command + "\n")
    f.close()
    shell(command)
    filetosend = ""

if snakemake.params.filesender_libsize:
    command = "tar -czvf " + snakemake.output.libsize + " " + snakemake.params.res_libsize + "* >> " + log_filename + " 2>&1"
    f = open(log_filename, 'at')
    f.write("## COMMAND: " + command + "\n")
    f.close()
    shell(command)
    filetosend = " " + snakemake.output.libsize
else:
    command = "touch " + snakemake.output.libsize + " >> " + log_filename + " 2>&1"
    f = open(log_filename, 'at')
    f.write("## COMMAND: " + command + "\n")
    f.close()
    shell(command)
    filetosend = ""

command = "python3 wrappers/filesender/filesender.py -m \"" + snakemake.params.message + "\" -s \"" + snakemake.params.subject + "\" -r \"" + snakemake.params.recipient + "\" " + filetosend + " >> " + log_filename + " 2>&1" #print command without credentials
f = open(log_filename, 'at')
f.write("## COMMAND: " + command + "\n")
f.close()
shell("python3 wrappers/filesender/filesender.py -m \"" + snakemake.params.message + "\" -s \"" + snakemake.params.subject + "\" -u " + username + " -a " + apikey + " -r \"" + snakemake.params.recipient + "\" " + filetosend + " >> " + log_filename + " 2>&1")
