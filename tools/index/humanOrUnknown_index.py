"""
Function:     Generate human (or unknown) sequence index
"""

import os
from optparse import OptionParser

def scan_genome_path(file_path):
    """
    Args:
        file_path: str; a string denotes the genomes' storage path            
    """
    # obtain all of the absolute paths of microbial genome files
    file_path = file_path.rstrip('/') + '/'                                # in case no '/' is at the end of the input file_path
    files_abspath = []
    files = os.listdir(file_path)
    for file_ in files:
        if file_.endswith(".fna") or file_.endswith(".fa") or file_.endswith(".fasta"):
            file_ = file_path + file_
            files_abspath.append(file_)
    # scan all of the microbial genome files and build a index of them
    index_path = os.path.abspath(file_path + "../") + "/genome_index/"    # make a new folder in the parent directory of microbial
                                                                        # genome storage path to store the microbial genome index file
    try:
        os.makedirs(index_path)
    except:
        pass                                                                    
    genome_index_abspath = index_path + "noisy_data_index"
    with open(genome_index_abspath, 'w') as f_w:
        for file_ in files_abspath:
            genome_length = 0
            with open(file_) as f_r:
                annotation = f_r.readline().strip()
                name = annotation.lstrip('>').split(' ', 1)[0]
                if ',' in name:
                    name = name.split(',')[0]
                for line in f_r:
                    line = line.strip()
                    genome_length += len(line)
            genome_info = name + '\t' + str(genome_length) + '\t' + file_
            f_w.write(genome_info)
            f_w.write('\n')
    print "Index file saved in %s." % genome_index_abspath

def main():
    usageReminder = "usage: python humanOrUnknown_index.py [-d] ref_seq_folder"
    parser = OptionParser(usage = usageReminder)
    parser.add_option("-d", action = "store", type = "string", dest = "ref_seq_folder", \
                      help = "The folder path of the human (or unknown) reference sequences to be used.")
    (opts, args) = parser.parse_args()
    genome_path = opts.ref_seq_folder
    scan_genome_path(genome_path)

if __name__ == '__main__':
    main()
