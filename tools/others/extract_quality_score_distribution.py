"""
Function:    Extract the quality score distribution of a fastq sequencing file
"""

import os
import sys
import numpy as np
from optparse import OptionParser

RANGES = {
    "ILLUMINA_ASCII_BASE=33": (33, 75),
    "ILLUMINA_ASCII_BASE=64": (64, 106)
}

def get_qual_range(qual_str):
    """
    Args:
        qual_str:         str; the ascii-encoded quality score
    """
    min_val, max_val = sys.maxint, -sys.maxint
    count = 0
    for char in qual_str:
        count += 1
        char = ord(char)
        if char < min_val:
            min_val = char
        if char > max_val:
            max_val = char
    return (min_val, max_val, count)

def get_encoding(guess_min, guess_max, ranges = RANGES):
    """
    Args:
        guess_min:         int; the minimum quality score of the input fastq file
        guess_max:        int; the maximum quality score of the input fastq file
        ranges:            dict; the quality score encoding methods  
    """
    valid_encoding = ""
    for encoding, (encode_min, encode_max) in ranges.items():
        if guess_min >= encode_min and guess_max <= encode_max:
            valid_encoding = encoding
    return valid_encoding

def check_ascii_coding(file_name):
    """
    Args:
        file_name:         str; fastq file
    """
    # count the number of total lines   
    line_count = 0
    with open(file_name) as f_r:
        for line in f_r:
            line_count += 1
    # check the quality score encoding method
    line_count = 0
    read_len_set = set()
    global_min, global_max = sys.maxint, -sys.maxint
    with open(file_name) as f_r:
         for line in f_r:
            line_count += 1
            line = line.rstrip()
            if line_count % 4 == 0:
                cur_min, cur_max, read_len = get_qual_range(line)
                read_len_set.add(read_len)
                if cur_min < global_min or cur_max > global_max:
                    global_min, global_max = min(cur_min, global_min), max(cur_max, global_max)
    encoding = get_encoding(global_min, global_max)
    max_read_len = max(read_len_set)
    if len(encoding) == 0:
        print "Invalid fastq file! No suitable quality score encoding method was detected!"
        os._exit(1)
    return encoding, max_read_len

def extract_qual_distri(file_name, encoding, max_read_len):
    """
    Args:
        file_name:         str; fastq file
        encoding:        str; the quality score encoding method, which can be searched in the global
                        varaiable "RANGES"
        max_read_len:    int; the maximum read length of the input fastq file 
    """
    # calculate the quality score distribution of each position of a read
    line_count = 0
    encode_min, encode_max = RANGES[encoding]
    encode_range = encode_max - encode_min + 1
    qual_score_tot_distri = np.zeros((max_read_len, encode_range))
    with open(file_name) as f_r:
        for line in f_r:
            line_count += 1
            line = line.rstrip()
            if line_count % 4 == 0:
                for i in xrange(len(line)):
                    score = ord(line[i]) - encode_min
                    qual_score_tot_distri[i][score] += 1
    # adjust the zero-count position with a small value to avoid zero probability
    filling_value = 0.01
    qual_score_tot_distri[qual_score_tot_distri == 0] = filling_value
    # normalize the distribution
    line_sum = qual_score_tot_distri.sum(axis = 1)
    line_sum = np.reshape(line_sum, (max_read_len, 1))
    line_sum = np.tile(line_sum, [1, encode_range])
    qual_score_tot_distri = qual_score_tot_distri / line_sum
    # calculate the cumulative probability distribution
    qual_score_tot_distri = qual_score_tot_distri.cumsum(axis = 1)
    return qual_score_tot_distri

def write_qual_distri(qual_score_distri, encoding, file_path):
    """
    Args:
        qual_score_dis:    np.array; the quality score distribution calculated by function "extract_qual_distri"
        encoding:        str; the encoding method of quality score
        file_path:        str; the absolute path to save the quality score distribution file
    """
    file_path = file_path.rstrip('/') + '/'
    file_name_write = file_path + "quality_score.txt"
    max_read_len, encode_range = qual_score_distri.shape
    with open(file_name_write, 'w') as f_w:
        f_w.write("#quality score probability distribution")
        f_w.write('\n')
        f_w.write("#%s" % encoding)
        f_w.write('\n')
        f_w.write("#score_range=0-%s" % str(encode_range - 1))
        f_w.write('\n')
        f_w.write("#max_read_length=%s" % str(max_read_len))
        f_w.write('\n')
        for i in xrange(max_read_len):
            qual_score_distri_read = list(qual_score_distri[i])
            qual_score_distri_read = ':'.join([str(s) for s in qual_score_distri_read])
            f_w.write("pos_%s\t%s" % (str(i), qual_score_distri_read))
            f_w.write('\n')
    print "Quality score distribution file saved in %s." % file_name_write

def main():
    usageReminder = ("usage: python extract_quality_score_distribution.py [-i] input_real_seq_data "
                     "[-o] output_qual_distri_dir")
    parser = OptionParser(usage = usageReminder)
    parser.add_option("-i", action = "store", type = "string", dest = "input_real_seq_data", \
                      help = "Input sequencing data for quality score extraction.")
    parser.add_option("-o", action = "store", type = "string", dest = "output_qual_distri_dir", \
                      help = "Output folder path for saving the extracted quality score distribution.")
    (opts, args) = parser.parse_args()
    file_name = opts.input_real_seq_data
    file_path = opts.output_qual_distri_dir
#    file_name = "/data/ssliu/Simulation/SimExperimentPipe/simMicrobeEnv/genNonuniReads/real_seq_data/DLF001.pair.1.fq"
#    file_path = "/data/ssliu/Simulation/SimExperimentPipe/simMicrobeEnv/genNonuniReads/real_seq_data/"

    encoding, max_read_len = check_ascii_coding(file_name)
    qual_score_tot_distri = extract_qual_distri(file_name, encoding, max_read_len)
    write_qual_distri(qual_score_tot_distri, encoding, file_path)

if __name__ == '__main__':
    main()
