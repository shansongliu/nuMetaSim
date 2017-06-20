"""
Function:     calculate the relationship between GC content and coverage bias 
"""

import os
import random
from optparse import OptionParser

def list_genome(file_path):
    file_path = file_path.rstrip('/') + '/'
    file_abs_paths = []
    file_names = os.listdir(file_path)
    for file_ in file_names:
        if file_.endswith(".fna"):
            file_name = file_path + file_
            file_abs_paths.append(file_name)
    return file_abs_paths

def read_genome(file_name):
    annotation = ""
    genome = []
    start_end = []
    line_count = 0
    with open(file_name) as f_r:
        annotation = f_r.readline().strip()
        start, end = 0, 0
        for line in f_r:
            line_count += 1
            line = line.strip()
            line_len = len(line)
            if line_count == 1:
                start, end = 0, start + line_len
            else:
                start, end = end, end + line_len
            genome.extend(line)
            start_end.append((start, end))
    return annotation, genome, start_end

def write_shuffled_genome(annotation, shuffled_genome, start_end, file_name_write):
    with open(file_name_write, 'w') as f_w:
        annot = annotation + "|shuffled"
        f_w.write(annot)
        f_w.write('\n')
        for pos in start_end:
            write_line = shuffled_genome[pos[0]:pos[1]]
            write_line = write_line.strip()
            f_w.write(write_line)
            f_w.write('\n')
    print ("Shuffled genome saved in %s" % file_name_write)

def main():
    usageReminder = ("usage: python genome_shuffle.py [-d] ref_genome_folder")
    parser = OptionParser(usage = usageReminder)
    parser.add_option("-d", action = "store", type = "string", dest = "ref_genome_folder", \
                      help = "The folder path of all the genomes to be shuffled.")
    (opts, args) = parser.parse_args()

    genome_path = opts.ref_genome_folder
    genome_paths = list_genome(genome_path)
    for geo_path in genome_paths:
        annotation, genome, start_end = read_genome(geo_path)
        random.shuffle(genome)
        genome = ''.join(genome)
        shuffled_genome_path = genome_path.rstrip('/') + '/' + geo_path.split('/')[-1].rstrip(".fna") + "_shuffled.fna" 
        write_shuffled_genome(annotation, genome, start_end, shuffled_genome_path)

if __name__ == '__main__':
    main()
