######################################
# wrapper for rule: biobloom
######################################
import subprocess
import pandas as pd

from snakemake.shell import shell
shell.executable("/bin/bash")
log_filename = str(snakemake.log)
ref_dir = str(snakemake.params.ref_dir)
ref_list = str(snakemake.params.ref_list)
ref_table = ref_dir + "/reference_info/biobloom_reference_sequia.txt"

f = open(log_filename, 'a+')
f.write("\n##\n## RULE: biobloom \n##\n")
f.close()

version = str(subprocess.Popen("conda list ", shell=True, stdout=subprocess.PIPE).communicate()[0], 'utf-8')
f = open(log_filename, 'at')
f.write("## CONDA: "+version+"\n")
f.close()

#Minimum aligned reads to run Biobloom (%):

# command = "cat " + snakemake.input.flagstats + " | grep -P \"^[0-9]+ \+ [0-9]+ mapped\" | sed \"s/[0-9]\+ + [0-9]\+ mapped (\([^%]\+\)% .\+/\\1/\" "
# with open(log_filename, 'at') as f:
#     f.write("## COMMAND: " + command + "\n")
#
# mapped_reads = str(subprocess.Popen(command,shell=True,stdout=subprocess.PIPE).communicate()[0], 'utf-8')
# f = open(log_filename, 'at')
# f.write("## Number of mapped reads: "+str(mapped_reads)+"\n")
# f.close()
#
# if float(mapped_reads) > snakemake.params.max_mapped_reads_to_run:

# set up contamination filters
biobloom_ref = pd.read_csv(ref_table, sep=";")
biobloom_ref["path"] = ref_dir + "/" + biobloom_ref["path"]
ref_list = ref_list.split(";")
filters_list = list()
for r in ref_list:
    if r not in biobloom_ref["species"].unique():
        with open(log_filename, 'at') as l:
            l.write("## WARNING: " + r + " is not a proper reference name. Reference must be in the form 'H.sapiens' \n")

    else:
        path = biobloom_ref[biobloom_ref["species"].str.contains(r)]["path"].values[0]
        filters_list.append(path)



# human_38 = ref_dir + "/homo_sapiens/GRCh38-p10/index/BioBloomTools/human_38.bf"
# mouse    = ref_dir + "/mus_musculus/GRCm38.p6-93/index/BioBloomTools/mouse.bf"
# #rat      = ref_dir + "/rattus_norvegicus/Rnor_6.0-91/index/BioBloomTools/rat.bf"
# yeast    = ref_dir + "/saccharomyces_cerevisiae/R64-1-1.100/index/BioBloomTools/yeast.bf"
# fruit_fly= ref_dir + "/drosophila_melanogaster/BDGP6-99/index/BioBloomTools/fruit_fly.bf"
# #dog      = ref_dir + "/canis_familiaris/CanFam3.1-101/index/BioBloomTools/dog.bf"
# athaliana= ref_dir + "/arabidopsis/TAIR10-31/index/BioBloomTools/A.thaliana.bf"
# #brassica = ref_dir + "/brassica_napus/Bra_napus_v2/index/BioBloomTools/B.napus.bf"
# c_elegans= ref_dir + "/c_elegans/WBcel235-102/index/BioBloomTools/C.elegans.bf"

# filters_list = list()
# if snakemake.params.filters == "all":
#     filters_list.append(human_38)
#     filters_list.append(mouse)
#     #filters_list.append(rat)
#     filters_list.append(yeast)
#     filters_list.append(fruit_fly)
#     #filters_list.append(dog)
#     filters_list.append(athaliana)
#     #filters_list.append(brassica)
#     filters_list.append(c_elegans)
# else:
#     for f in filters.split(';'):
#         if f == 'human':
#             filters_list.append(human_38)
#         elif f == 'mouse':
#             filters_list.append(mouse)
#         # elif f == 'rat':
#         #     filters_list.append(rat)
#         elif f == 'yeast':
#             filters_list.append(yeast)
#         elif f == 'fly':
#             filters_list.append(fruit_fly)
#         # elif f == 'dog':
#         #     filters_list.append(dog)
#         elif f == 'arabidopsis':
#             filters_list.append(athaliana)
#         # elif f == 'brassica':
#         #     filters_list.append(brassica)
#         elif f == 'c_elegans':
#             filters_list.append(c_elegans)
#         else:
#             with open(log_filename, 'at') as l:
#                 l.write("## WARNING: non existing species used for contamination check: "+f+"\n")

command = 'mkdir -p $(dirname '+snakemake.output.table+') >> '+log_filename+" 2>&1"
f = open(log_filename, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)


if snakemake.params.paired == "SE":
    with open(log_filename, 'at') as f:
        f.write("## INFO: I noticed Single-end data\n")

    command = "(time biobloomcategorizer -p " + snakemake.params.prefix + \
              " -t " + str(snakemake.threads) + \
              " -f '" + " ".join(filters_list) + "'" + \
              " <(zcat " + snakemake.input.r1 + ")" + \
              " ) >> " + log_filename + " 2>&1 "
elif snakemake.params.paired == "PE":
    with open(log_filename, 'at') as f:
        f.write("## INFO: I noticed Paired-end data\n")

    command = "(time biobloomcategorizer -p " + snakemake.params.prefix + \
              " -t " + str(snakemake.threads) + \
              " -e -f '" + " ".join(filters_list) + "'" + \
              " <(zcat " + snakemake.input.r1 + ")" + \
              " <(zcat " + snakemake.input.r2 + ")" + \
              " ) >> " + log_filename + " 2>&1 "
else:
    with open(log_filename, 'at') as f:
        f.write("## WARNING: something went wrong, there is no input in snakemake.input object\n")

f = open(log_filename, 'at')
f.write("## COMMAND: " + command + "\n")
f.close()
shell(command)

# else:
#     command = "cat \"Not enough aligned reads ("+ mapped_reads +"%). Not neccesary to run Biobloom.\n\" > " + snakemake.output.table
#     f = open(log_filename, 'at')
#     f.write("## COMMAND: " + command + "\n")
#     f.close()
#     shell(command)

