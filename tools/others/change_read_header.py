"""
Function:    Change the description line of each read
"""

import os
from optparse import OptionParser

usageReminder = ("usage: python [-i] input_file_name [-o] output_file_name [-l] fa_fq_flag")
parser = OptionParser(usage = usageReminder)
parser.add_option("-i", action = "store", type = "string", dest = "input_file_name", \
                  help = "Name of the input file (or absolute path of the file).")
parser.add_option("-o", action = "store", type = "string", dest = "output_file_name", \
                  help = "Name of the output file (or absolute path of the file).")
parser.add_option("-l", action = "store", type = "int", dest = "fa_fq_flag", \
                  help = "Fasta or fastq flag, denoting the input file format.")
(opts, args) = parser.parse_args()

input_file = opts.input_file_name
output_file = opts.output_file_name
fa_fq_flag = opts.fa_fq_flag

if fa_fq_flag == None:
    print "Fasta or fastq flag is needed. Program will exit."
    os._exit(1)
elif fa_fq_flag == 0:
    if not (input_file.endswith(".fa") or input_file.endswith(".fasta")):
        print "The fasta or fastq flag is set as 0, the input file's suffix should be \".fa\" or \".fasta\". Program with exit."
        os._exit(1)
elif fa_fq_flag == 1:
    if not (input_file.endswith(".fq") or input_file.endswith(".fastq")):
        print "The fasta or fastq flag is set as 0, the input file's suffix should be \".fq\" or \".fastq\". Program with exit."
        os._exit(1) 

count = 0
read_name_tmp = ""
with open(output_file, 'w') as f_w:
    with open(input_file) as f_r:
        for line in f_r:
            if line.startswith('>') or line.startswith("@>"):
                count += 1
                anno = line.split()
                if fa_fq_flag == 0:
                    read_name = '>' + anno[0].lstrip('>')
                elif fa_fq_flag == 1:
                    read_name = '@>' + anno[0].lstrip('@>')                    
                if read_name_tmp != read_name:
                    count = 0
                    read_name_tmp = read_name
                read_name = read_name + "_%s" % str(count)
                anno = [read_name, ' '.join(anno[1:])]
                anno = ' '.join(anno)
                f_w.write(anno)
                f_w.write('\n')
            else:
                f_w.write(line)
