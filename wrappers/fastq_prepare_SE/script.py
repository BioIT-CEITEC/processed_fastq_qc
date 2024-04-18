######################################
# wrapper for rule: fastq_prepare_PE
######################################
import os
import subprocess
import sys
from snakemake.shell import shell
import json

log_filename = str(snakemake.log)

f = open(log_filename, 'wt')
f.write("\n##\n## RULE: fastq_prepare_PE \n##\n")
f.close()

shell.executable("/bin/bash")

command = "mkdir -p " + os.path.dirname(snakemake.output[0])
f = open(log_filename, 'at')
f.write("## CREATE_OUTPUT_DIR: " + command + "\n")
f.close()
shell(command)


def load_and_configure_UMI(wf_config_path,config):
    with open(wf_config_path, 'r') as file:
        wf_config = json.load(file)

    # Extract the primary GUI parameters from the config
    primary_gui_params = wf_config['gui_params']['primary']

    # Initialize the dictionary to store UMI settings
    umi_settings = {}

    # The list of parameters we're interested in setting based on UMI type
    param_keys = [
        "UMI_write_to","UMI_R1_start", "UMI_R1_end", "insert_R1_start",
        "UMI_R2_start", "UMI_R2_end", "insert_R2_start"
    ]

    # Extract values for each UMI type under each parameter key
    for umi_type in primary_gui_params['UMI']['list'].keys():
        settings = {}
        for param_key in param_keys:
            if 'conditions' in primary_gui_params[param_key] and \
                    'value' in primary_gui_params[param_key]['conditions'] and \
                    'UMI' in primary_gui_params[param_key]['conditions']['value'] and \
                    umi_type in primary_gui_params[param_key]['conditions']['value']['UMI']:
                settings[param_key] = primary_gui_params[param_key]['conditions']['value']['UMI'][umi_type]

        # Store settings for this UMI type if we have any settings defined
        if settings:
            umi_settings[umi_type] = settings

    # Get the UMI type from config
    umi_type = config.get("UMI")

    # Set the config values based on the UMI type
    if umi_type in umi_settings:
        for key, value in umi_settings[umi_type].items():
            config[key] = value

    return config



def replace_last_occurrence(s, old, new):
    # Reverse the string (so the last occurrence becomes the first)
    reversed_s = s[::-1]

    # Replace the first occurrence of the substring in the reversed string
    reversed_s = reversed_s.replace(old[::-1], new[::-1], 1)

    # Reverse the string back to its original orientation
    s = reversed_s[::-1]

    return s

def get_header_separator_character(in_filename):
    command = "zcat " + in_filename + " | head -n 1"

    # Execute the command
    output = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE).communicate()[0]
    first_line = str(output, 'utf-8')

    # Check if the first line contains " " or "/"
    if " " in first_line:
        sep = " "
    elif "/" in first_line:
        sep = "/"

    return sep


if os.stat(snakemake.input.in_filename).st_size != 0:
    sample = snakemake.wildcards.sample
    in_filename = snakemake.input.in_filename
    in_filename_R2 = replace_last_occurrence(in_filename, "_R1", "_R2")
    config = load_and_configure_UMI("workflow.config.json",snakemake.params.config)
    umi = config["UMI"]
    print("processing:" + sample + " with UMI: " + umi)

    sep = get_header_separator_character(in_filename)

    unzip_out_R1 = snakemake.output.fastq[:-3]

    R1_UMI_start = config["UMI_R1_start"]
    R1_UMI_end = config["UMI_R1_end"]
    R1_ins_start = config["insert_R1_start"]
    R2_UMI_start = config["UMI_R2_start"]
    R2_UMI_end = config["UMI_R2_end"]

    if R1_ins_start < 1:
        R1_ins_start = 1

    if R1_UMI_start > R1_UMI_end:
        R1_UMI_start = 0
    if R2_UMI_start > R2_UMI_end:
        R2_UMI_start = 0

    if umi != "no_umi" and (config["UMI_from_R3_file"] or (R1_UMI_start != 0 or R2_UMI_start != 0)):
        if config["UMI_from_R3_file"]:
            in_filename_UMI = replace_last_occurrence(in_filename, "_R1", "_UMI")
            if not os.path.isfile(in_filename_UMI):
                sys.exit("UMI from extra fastq specified, but it doesn't exist.")

            if config["UMI_write_to"] == "sep_file":
                umi_file = os.path.dirname(snakemake.output.R1) + "/" + sample + ".UMI.fastq"
                command = "gunzip -c " + in_filename_UMI + " > " + umi_file + "; " +\
                          "mv -T " + in_filename + " " + snakemake.output.fastq

            elif config["UMI_write_to"] == "sep_file_gz":
                umi_file = os.path.dirname(snakemake.output.R1) + "/" + sample + ".UMI.fastq.gz"
                command = "mv " + in_filename_UMI + " " + umi_file + "; " +\
                          "mv -T " + in_filename + " " + snakemake.output.fastq
            else:
                command = "(paste <(zcat " + in_filename + ") <(zcat " + in_filename_UMI + ") |" + \
                          " awk '{{ if(NR%4==1) {{split($1,head_R1,\"" + sep + "\")" + \
                          " else if(NR%4==2) {{print head_R1[1] \"_\" $2 \" \" head_R1[2] \"\\n\" substr($1,7) > out1" + \
                          "}} else if(NR%4==0) {{print substr($1,7) > out1" + \
                          "}} else {{print $1 > out1}} }}' FS='\\t' out1=" + unzip_out_R1 + \
                          " && gzip -f " + unzip_out_R1 + ") 2>> " + log_filename
        else:
            if R1_UMI_start != 0:
                sub_R1_UMI = "substr($1, "+ str(R1_UMI_start) +","+str(R1_UMI_end)+")"
            else:
                sub_R1_UMI = ""

            if R2_UMI_start != 0:
                if not os.path.isfile(in_filename_R2):
                    sys.exit("UMI from R2 fastq specified, but R2 fastq file does not exist.")
                sub_R2_UMI = "substr($1, " + str(R2_UMI_start) + "," + str(R2_UMI_end) + ") "
                in_file_R2_text = "<(zcat " + in_filename_R2 + ")"
            else:
                sub_R2_UMI = ""
                in_file_R2_text = ""

            if config["UMI_write_to"] == "fastq_header":
                command = "(paste <(zcat " + in_filename + ") " + in_file_R2_text + "|" + \
                          " awk '{{ if(NR%4==1) {{split($1,head_R1,\"" + sep + "\")}}" + \
                          " else if(NR%4==2) {{umi=" + sub_R1_UMI + sub_R2_UMI + "; print head_R1[1] \"_\" umi \" \" head_R1[2] \"\\n\" substr($1," + str(R1_ins_start) + ") > out1" + \
                          "}} else if(NR%4==0) {{print substr($1," + str(R1_ins_start) + ") > out1" + \
                          "}} else {{print $1 > out1}} }}' FS='\\t' out1=" + unzip_out_R1 +  \
                          " && gzip -f " + unzip_out_R1 + ") 2>> " + log_filename
            elif config["UMI_write_to"] == "sep_file":
                umi_file = os.path.dirname(snakemake.output.R1) + "/" + sample + ".UMI.fastq"
                command = "(paste <(zcat " + in_filename + ") " + in_file_R2_text + "|" + \
                          " awk '{{ if(NR%4==1) {{split($1,head_R1,\"" + sep + "\")}}" + \
                          " else if(NR%4==2) {{umi=" + sub_R1_UMI + sub_R2_UMI + "; print head_R1[1] \"" + sep + "\" head_R1[2] \"\\n\" substr($1," + str(R1_ins_start) + ") > out1;" + \
                          " print head_R1[1] \"\\n\" umi > umi_out}} else if(NR%4==0) {{print substr($1," + str(R1_ins_start) + ") > out1;" + \
                          " print " + sub_R1_UMI + sub_R2_UMI + " > umi_out}} else {{print $1 > out1; print $1 > umi_out}} }}' FS='\\t' out1=" + unzip_out_R1 + " umi_out=" + umi_file + \
                          " && gzip -f " + unzip_out_R1 + ") 2>> " + log_filename
            else:
                umi_file = os.path.dirname(snakemake.output.R1) + "/" + sample + ".UMI.fastq"
                command = "(paste <(zcat " + in_filename + ") " + in_file_R2_text + "|" + \
                          " awk '{{ if(NR%4==1) {{split($1,head_R1,\"" + sep + "\")}}" + \
                          " else if(NR%4==2) {{umi=" + sub_R1_UMI + sub_R2_UMI + "; print head_R1[1] \"" + sep + "\" head_R1[2] \"\\n\" substr($1," + str(R1_ins_start) + ") > out1;" + \
                          " print head_R1[1] \"\\n\" umi > umi_out}} else if(NR%4==0) {{print substr($1," + str(R1_ins_start) + ") > out1;" + \
                          " print " + sub_R1_UMI + sub_R2_UMI + " > umi_out}} else {{print $1 > out1; print $1 > umi_out}} }}' FS='\\t' out1=" + unzip_out_R1 + " umi_out=" + umi_file + \
                          " && gzip -f " + unzip_out_R1 + " " + umi_file + ") 2>> " + log_filename



    else:
        command = "mv -T " + in_filename + " " + snakemake.output.fastq


else:
    command = "touch " + snakemake.output.fastq

f = open(log_filename, 'at')
f.write("## UMI COMMAND: " + command + "\n")
f.close()
shell(command)

    # elif umi == "IDT":
    #     out_R1 = gzip.open(snakemake.output.R1, 'wt')
    #     out_R2 = gzip.open(snakemake.output.R2, 'wt')
    #
    #     in_filename_UMI = in_filename_R2
    #     in_filename_R2 = re.sub("_R2_", "_R3_", in_filename_R2)
    #
    #     with gzip.open(in_filename, 'rt') as R1, gzip.open(in_filename_R2, 'rt') as R2, gzip.open(in_filename_UMI,
    #                                                                                               'rt') as UMI:
    #         i = 0
    #         for R1_line, R2_line, UMI_line in zip(R1, R2, UMI):
    #             i += 1
    #             if i % 4 == 1:
    #                 header_R1 = R1_line.strip()
    #                 header_R2 = R2_line.strip()
    #             elif i % 4 == 2:
    #                 out_R1.write(header_R1.split(" ")[0] + "_" + UMI_line.strip() + " " + header_R1.split(" ")[
    #                     1] + "\n" + R1_line)
    #                 out_R2.write(header_R2.split(" ")[0] + "_" + UMI_line.strip() + " " + header_R2.split(" ")[
    #                     1] + "\n" + R2_line)
    #             elif i % 4 == 0:
    #                 out_R1.write(R1_line)
    #                 out_R2.write(R2_line)
    #             else:
    #                 out_R1.write(R1_line)
    #                 out_R2.write(R2_line)
    #
    # elif umi == "BRONCO":
    #     out_R1 = gzip.open(snakemake.output.R1, 'wt')
    #     out_R2 = gzip.open(snakemake.output.R2, 'wt')
    #
    #     in_filename_UMI = in_filename_R2
    #     in_filename_R2 = re.sub("_R2_", "_R3_", in_filename_R2)
    #
    #     with gzip.open(in_filename, 'rt') as R1, gzip.open(in_filename_R2, 'rt') as R2, gzip.open(in_filename_UMI,
    #                                                                                               'rt') as UMI:
    #         i = 0
    #         for R1_line, R2_line, UMI_line in zip(R1, R2, UMI):
    #             i += 1
    #             if i % 4 == 1:
    #                 header_R1 = R1_line.strip()
    #                 header_R2 = R2_line.strip()
    #             elif i % 4 == 2:
    #                 out_R1.write(header_R1.split(" ")[0] + "_" + UMI_line.strip() + " " + header_R1.split(" ")[
    #                     1] + "\n" + R1_line)
    #                 out_R2.write(header_R2.split(" ")[0] + "_" + UMI_line.strip() + " " + header_R2.split(" ")[
    #                     1] + "\n" + R2_line)
    #             elif i % 4 == 0:
    #                 out_R1.write(R1_line)
    #                 out_R2.write(R2_line)
    #             else:
    #                 out_R1.write(R1_line)
    #                 out_R2.write(R2_line)
    #
    # elif umi == "LYNX":
    #     umi_file = os.path.dirname(snakemake.output.R2) + "/" + sample + ".UMI.fastq"
    #     umi_file_in = in_filename_R2
    #     in_filename_R2 = re.sub("_R2_", "_R3_", in_filename_R2)
    #
    #     command = "gunzip -c " + umi_file_in + " > " + umi_file
    #     f = open(log_filename, 'at')
    #     f.write("## UMI COMMAND: " + command + "\n")
    #     f.close()
    #     shell(command)
    #
    #     # copy R1
    #     command = "mv -T " + in_filename + " " + snakemake.output.R1
    #
    #     f = open(log_filename, 'at')
    #     f.write("## COMMAND: " + command + "\n")
    #     f.close()
    #     shell(command)
    #
    #     # copy R2
    #     command = "mv -T " + in_filename_R2 + " " + snakemake.output.R2
    #
    #     f = open(log_filename, 'at')
    #     f.write("## COMMAND: " + command + "\n")
    #     f.close()
    #     shell(command)
    #
    # elif umi == "Qiaseq":
    #     out_R1 = snakemake.output.R1[:-3]
    #     out_R2 = snakemake.output.R2[:-3]
    #
    #     command = "(paste <(zcat " + in_filename + ") <(zcat " + in_filename_R2 + ") |" + \
    #               " awk '{{ if(NR%4==1) {{split($1,head_R1,\" \"); split($2,head_R2,\" \")}}" + \
    #               " else if(NR%4==2) {{umi=substr($2,1,12); print head_R1[1] \"_\" umi \" \" head_R1[2] \"\\n\" $1 > out1;" + \
    #               " print head_R2[1] \"_\" umi \" \" head_R2[2] \"\\n\" substr($2,24) > out2}} else if(NR%4==0) {{print $1 > out1;" + \
    #               " print substr($2,24) > out2}} else {{print $1 > out1; print $2 > out2}} }}' FS='\\t' out1=" + out_R1 + " out2=" + out_R2 + \
    #               " && gzip -f " + out_R1 + " " + out_R2 + ") 2>> " + log_filename
    #     with open(log_filename, 'at') as f:
    #         f.write("## COMMAND: " + command + "\n")
    #     shell(command)

        # if "externally_sequenced_fake" in run_name:
        #    command = "cp -T "+in_filename+" "+ snakemake.output.R1
        # else:
        #    command = "mv -T "+in_filename+" "+ snakemake.output.R1

        # f = open(log_filename, 'at')
        # f.write("## COMMAND: "+command+"\n")
        # f.close()
        # shell(command)

        # umi_file = os.path.dirname(snakemake.output.R2)+ "/"+sample+".UMI.fastq"
        # shell("gunzip "+in_filename_R2)
        # in_filename_R2 = in_filename_R2.replace(".gz","")

        # f_R2 = open(snakemake.output.R2.replace(".gz",""), 'w')
        # f_umi = open(umi_file, 'w')

        # with(open(in_filename_R2,'r')) as f_in:
        #    for line_num,line in enumerate(f_in):
        #        if line_num % 2 == 1:
        #            f_umi.write(line[:12]+"\n")
        #            f_R2.write(line[24:])
        #        else:
        #            f_umi.write(line)
        #            f_R2.write(line)

        # f_umi.close()
        # f_R2.close()
        # shell("gzip " + snakemake.output.R2.replace(".gz",""))
        # shell("rm "+in_filename_R2)
#
#     elif umi == "CS_UMI_sep_file":
#         umi_file = os.path.dirname(snakemake.output.R2) + "/" + sample + ".UMI.fastq"
#         out_R1 = snakemake.output.R1[:-3]
#         out_R2 = snakemake.output.R2[:-3]
#
#         # Command to decompress the first line of the file and capture the output
#         command = "zcat " + in_filename + " | head -n 1"
#
#         # Execute the command
#         output = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE).communicate()[0]
#         first_line = str(output, 'utf-8')
#
#         # Check if the first line contains " " or "/"
#         if " " in first_line:
#             sep = " "
#         elif "/" in first_line:
#             sep = "/"
#
#         command = "(paste <(zcat " + in_filename + ") <(zcat " + in_filename_R2 + ") |" + \
#                   " awk '{{ if(NR%4==1) {{split($1,head_R1,\"" + sep + "\"); split($2,head_R2,\"" + sep + "\")}}" + \
#                   " else if(NR%4==2) {{umi=substr($1,1,3)substr($2,1,3); print head_R1[1] \"" + sep + "\" head_R1[2] \"\\n\" substr($1,7) > out1;" + \
#                   " print head_R2[1] \"" + sep + "\" head_R2[2] \"\\n\" substr($2,7) > out2; print head_R1[1] \"\\n\" umi > umi_out}} else if(NR%4==0) {{print substr($1,7) > out1;" + \
#                   " print substr($2,7) > out2; print substr($1,1,3)substr($2,1,3) > umi_out}} else {{print $1 > out1; print $2 > out2; print $1 > umi_out}} }}' FS='\\t' out1=" + out_R1 + " out2=" + out_R2 + " umi_out=" + umi_file + \
#                   " && gzip -f " + out_R1 + " " + out_R2 + ") 2>> " + log_filename
#         with open(log_filename, 'at') as f:
#             f.write("## COMMAND: " + command + "\n")
#         shell(command)
#
#         # umi_file = os.path.dirname(snakemake.output.R2) + "/" + sample + ".UMI.fastq"
#         # out_R1 = snakemake.output.R1[:-3]
#         # out_R2 = snakemake.output.R2[:-3]
#         #
#         # command = "(paste -d '' <(zcat " + in_filename + " | awk '{{ if(NR%4==2 || NR%4==0) {{print substr($0,1,3); print substr($0,7) > out}} else {{print $0; print $0 > out}} }}' out=" + out_R1 + ") <(zcat " + in_filename_R2 + " | awk '{{ if(NR%4==2 || NR%4==0) {{print substr($0,1,3); print substr($0,7) > out}} else {{print \"\"; print $0 > out}} }}' out=" + out_R2 + ") > " + umi_file + " && gzip -f " + out_R1 + " " + out_R2 + ") 2>> " + log_filename
#         # with open(log_filename, 'at') as f:
#         #     f.write("## COMMAND: " + command + "\n")
#         # shell(command)
#         #
#         # command = "(paste <(zcat " + in_filename + ") <(zcat " + in_filename_R2 + ") |" + \
#         #           " awk '{{ if(NR%4==1) {{split($1,head_R1,\"/\"); split($2,head_R2,\"/\")}}" + \
#         #           " else if(NR%4==2) {{umi=substr($1,1,3)substr($2,1,3); print head_R1[1] \"_\" umi \" \" head_R1[2] \"\\n\" substr($1,7) > out1;" + \
#         #           " print head_R2[1] \"_\" umi \" \" head_R2[2] \"\\n\" substr($2,7) > out2}} else if(NR%4==0) {{print substr($1,7) > out1;" + \
#         #           " print substr($2,7) > out2}} else {{print $1 > out1; print $2 > out2}} }}' FS='\\t' out1=" + out_R1 + " out2=" + out_R2 + \
#         #           " && gzip -f " + out_R1 + " " + out_R2 + ") 2>> " + log_filename
#         # with open(log_filename, 'at') as f:
#         #     f.write("## COMMAND: " + command + "\n")
#         # shell(command)
#
#
#     elif umi == "CS_UMI":
#         out_R1 = snakemake.output.R1[:-3]
#         out_R2 = snakemake.output.R2[:-3]
#
#         # Command to decompress the first line of the file and capture the output
#         command = "zcat " + in_filename + " | head -n 1"
#
#         # Execute the command
#         output = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE).communicate()[0]
#         first_line = str(output, 'utf-8')
#
#         # Check if the first line contains " " or "/"
#         if " " in first_line:
#             sep = " "
#         elif "/" in first_line:
#             sep = "/"
#
#         command = "(paste <(zcat " + in_filename + ") <(zcat " + in_filename_R2 + ") |" + \
#                   " awk '{{ if(NR%4==1) {{split($1,head_R1,\"" + sep + "\"); split($2,head_R2,\"" + sep + "\")}}" + \
#                   " else if(NR%4==2) {{umi=substr($1,1,3)substr($2,1,3); print head_R1[1] \"_\" umi \" \" head_R1[2] \"\\n\" substr($1,7) > out1;" + \
#                   " print head_R2[1] \"_\" umi \" \" head_R2[2] \"\\n\" substr($2,7) > out2}} else if(NR%4==0) {{print substr($1,7) > out1;" + \
#                   " print substr($2,7) > out2}} else {{print $1 > out1; print $2 > out2}} }}' FS='\\t' out1=" + out_R1 + " out2=" + out_R2 + \
#                   " && gzip -f " + out_R1 + " " + out_R2 + ") 2>> " + log_filename
#         with open(log_filename, 'at') as f:
#             f.write("## COMMAND: " + command + "\n")
#         shell(command)
#
#     elif umi == "TruSight_Oncology":
#         out_R1 = snakemake.output.R1[:-3]
#         out_R2 = snakemake.output.R2[:-3]
#
#         command = "(paste <(zcat " + in_filename + ") <(zcat " + in_filename_R2 + ") |" + \
#                   " awk '{{ if(NR%4==1) {{split($1,head_R1,\" \"); split($2,head_R2,\" \")}}" + \
#                   " else if(NR%4==2) {{umi=substr($1,1,6)substr($2,1,6); print head_R1[1] \"_\" umi \" \" head_R1[2] \"\\n\" substr($1,10) > out1;" + \
#                   " print head_R2[1] \"_\" umi \" \" head_R2[2] \"\\n\" substr($2,10) > out2}} else if(NR%4==0) {{print substr($1,10) > out1;" + \
#                   " print substr($2,10) > out2}} else {{print $1 > out1; print $2 > out2}} }}' FS='\\t' out1=" + out_R1 + " out2=" + out_R2 + \
#                   " && gzip -f " + out_R1 + " " + out_R2 + ") 2>> " + log_filename
#         with open(log_filename, 'at') as f:
#             f.write("## COMMAND: " + command + "\n")
#         shell(command)
#
#     else:
#         command = "mv -T " + in_filename + " " + snakemake.output.R1
#         f = open(log_filename, 'at')
#         f.write("## COMMAND: " + command + "\n")
#         f.close()
#         shell(command)
#
#         command = "mv -T " + in_filename_R2 + " " + snakemake.output.R2
#         f = open(log_filename, 'at')
#         f.write("## COMMAND: " + command + "\n")
#         f.close()
#         shell(command)
# else:
#     command = "touch " + snakemake.output.R1
#     f = open(log_filename, 'at')
#     f.write("## COMMAND: " + command + "\n")
#     f.close()
#     shell(command)
#
#     command = "touch " + snakemake.output.R2
#     f = open(log_filename, 'at')
#     f.write("## COMMAND: " + command + "\n")
#     f.close()
#     shell(command)
