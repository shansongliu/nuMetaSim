"""
Function:     Get the mapping files of the input sequencing reads to all reference genomes
"""

import os
import numpy as np
import multiprocessing
from optparse import OptionParser

def list_genome(file_path):
    """
    Args:
        file_path:         str; path of the input genomes                
    """
    file_path = file_path.rstrip('/') + '/'
    file_abs_paths = []
    file_names = os.listdir(file_path)
    for file_ in file_names:
        if file_.endswith(".fna"):
            file_name = file_path + file_
            file_abs_paths.append(file_name)
    return file_abs_paths

def split_genome(file_name, bin_size = 2000):
    """
    Args:
        file_name:         str; the genome file, it should only contains a single complete sequence
                        the file format is .fna (using .fa for short), which contains:
                            >U00096.3 Escherichia coli str. K-12 substr. MG1655, complete genome
                            ATCGAGTCAGTCAGTCGATCG...
                            CAGTCACGATAGTATCGGTCG...
                            GTCAGCAAGTCTCGATCGAGT...                
    """
    # read one genome
    annotation = ""
    genome = []
    with open(file_name) as f_r:
        annotation = f_r.readline().strip()        # e.g., >U00096.3 Escherichia coli str. K-12 substr. MG1655, complete genome
        acc_num = annotation.split()[0].lstrip('>')
        for line in f_r:
            genome.extend(line.strip())
    genome = ''.join(genome)

    # split the genome using bin_size
    genome_len = len(genome)
    n_bins = np.int(np.ceil(genome_len / np.float(bin_size)))
    file_path = '/'.join(file_name.split('/')[:-1]) + '/'
    file_name = file_name.split('/')[-1].rstrip(".fna") + "_split.fna"
    file_name = file_path + file_name
    with open(file_name, 'w') as f_w:
        for i in xrange(n_bins - 1):
            annot = '>' + str(i) + '-' + acc_num
            f_w.write(annot)
            f_w.write('\n')
            seq = genome[i * bin_size:(i + 1) * bin_size]
            f_w.write(seq)
            f_w.write('\n')
        annot = '>' + str(n_bins - 1) + '-' + acc_num
        f_w.write(annot)
        f_w.write('\n')
        seq = genome[(n_bins - 1) * bin_size:]
        f_w.write(seq)
        f_w.write('\n')
    return annotation

def list_split_genome(file_path):
    """
    Args:
        file_path:         str; path of the split genomes                
    """
    file_path = file_path.rstrip('/') + '/'
    file_abs_paths = []
    file_names = os.listdir(file_path)
    for file_ in file_names:
        if file_.endswith("_split.fna"):
            file_name = file_path + file_
            file_abs_paths.append(file_name)
    return file_abs_paths

def build_bowtie2_index(file_name):
    """
    Args:
        file_name:         str; the split genome file                
    """
    file_path = '/'.join(file_name.split('/')[:-1]) + '/'
    bowtie2_index = file_name.split('/')[-1].rstrip(".fna") + "_index"
    bowtie2_index = file_path + bowtie2_index
    os.system("bowtie2-build %s %s" % (file_name, bowtie2_index))
    return bowtie2_index

def map_with_bowtie2(file_name, reads_file):
    """
    Args:
        file_name:         str; the index prefix of each genome            
        reads_file:        str; input reads file name
    """
    file_path = '/'.join(file_name.split('/')[:-1]) + '/'
    split_genome_name_prefix = file_name.split('/')[-1].rstrip("_index")
    reads_file_name_prefix = reads_file.split('/')[-1].split('.')[0]
    sam_file_name = file_path + split_genome_name_prefix + '_' + reads_file_name_prefix + ".sam"
    os.system("bowtie2 -x %s -U %s -q -S %s --no-unal" % (file_name, reads_file, sam_file_name))

def list_sam_file(file_path):
    """
    Args:
        file_path:         str; path of the input sam files                
    """
    file_path = file_path.rstrip('/') + '/'
    file_abs_paths = []
    file_names = os.listdir(file_path)
    for file_ in file_names:
        if file_.endswith(".sam"):
            file_name = file_path + file_
            file_abs_paths.append(file_name)
    return file_abs_paths

def acc_microbe_name(file_path):
    """
    Args:
        file_path:         str; path of the input genomes, it should be the path for generating sam files                
    """    
    acc_micro_name_dict = {}
    file_path = file_path.rstrip('/') + '/'
    file_abs_paths = []
    file_names = os.listdir(file_path)
    for file_ in file_names:
        if file_.endswith(".fna") and (not file_.endswith("_split.fna")):
            file_ = file_path + file_
            f_r = open(file_)
            annotation = f_r.readline().strip()
            accession_num, microbe_name = annotation.lstrip('>').split(' ', 1)
            if ',' in microbe_name:
                microbe_name = microbe_name.split(',')[0]
            acc_micro_name_dict[accession_num] = microbe_name
            f_r.close()
    return acc_micro_name_dict

def read_count(file_name, acc_micro_name):
    """
    Args:
        file_name:         str; the .sam file to be counted
        acc_micro_name:    dict; {accession number: microbe_name}             
    """
    # count read in each genome bin
    read_count_dict = {}
    bin_name = ""
    with open(file_name) as f_r:
        for line in f_r:
            line = line.strip()
            if line.startswith('@'):
                if line.startswith("@SQ"):
                    bin_name = line.split()[1].split(':')[-1]
                    bin_num = np.int(bin_name.split('-')[0])
                    read_count_dict[bin_num] = 0
                else:
                    continue
            else:
                bin_name = line.split()[2]
                bin_num = np.int(bin_name.split('-')[0])
                read_count_dict[bin_num] += 1
    # normalize each genome bin using total count
    # to avoid 0 count genome bin, we add 0.01 to this kind of bin
    accession_num = bin_name.split('-')[-1]
    pseudo_count = 0.01
    total_mapped_count = 0
    for bin_num in read_count_dict:
        if read_count_dict[bin_num] == 0:
            read_count_dict[bin_num] = pseudo_count
            total_mapped_count += pseudo_count
        else:
            total_mapped_count += read_count_dict[bin_num]
    bin_propors = []
    for bin_num in read_count_dict:
        bin_propor = read_count_dict[bin_num] / np.float(total_mapped_count)
        bin_propors.append(str(bin_propor))
    bin_propors = ':'.join(bin_propors)
    return (acc_micro_name[accession_num], bin_propors)

def write_genome_distribution(file_path, geo_distri):
    """
    Args:
        file_path:         str; the absolute folder path to save the genome distribution file 
        geo_distri:        list; a list contains each calculated genome's distribution by counting reads 
                        which were mapped to bins in the genome            
    """
    file_path = file_path.rstrip('/') + '/'
    file_name = file_path + "distribution.txt"
    with open(file_name, 'w') as f_w:
        for term in geo_distri:
            microbe_name, distribution = term[0], term[1]
            write_line = microbe_name + '\t' + distribution
            f_w.write(write_line)
            f_w.write('\n')
    print "Genome distribution file saved in %s." % file_name

def main():
    usageReminder = ("usage: python cal_genome_distribution.py [--fq_path] fastq_data_path\n"
                     "[--geo_path] ref_genome_folder [-o] output_path [-b] bin_size")    
    parser = OptionParser(usage = usageReminder)
    parser.add_option("--fq_path", action = "store", type = "string", dest = "fastq_data_path", \
                      help = "The path of a real sequencing data sample. Absolute path is recommended.")
    parser.add_option("--geo_path", action = "store", type = "string", dest = "genome_path", \
                      help = "Input folder path for all the reference genomes.")
    parser.add_option("-o", action = "store", type = "string", dest = "output_path", \
                      help = "Output folder path for saving the genome distribution file.")
    parser.add_option("-b", action = "store", type = "int", dest = "bin_size", \
                      help = "The window size for mapping real sequencing reads. DeFault is 2000.")
    (opts, args) = parser.parse_args()
    if opts.bin_size == None:
        opts.bin_size = 2000

    real_data_path = opts.fastq_data_path
    genome_path = opts.genome_path
    genome_distribution_path = opts.output_path
    bin_size = opts.bin_size

    multiprocess_core = multiprocessing.cpu_count() / 2    # defalut processing core number (half of the total cpu cores)

    # split genome
    pool = multiprocessing.Pool(processes = multiprocess_core)
    genome_abs_paths = list_genome(genome_path)
    for i in xrange(len(genome_abs_paths)):
        pool.apply_async(split_genome, (genome_abs_paths[i], bin_size))
    pool.close()
    pool.join()

    # build bowtie2 index
    pool = multiprocessing.Pool(processes = multiprocess_core) 
    result = []
    split_genome_abs_paths = list_split_genome(genome_path)
    for i in xrange(len(split_genome_abs_paths)):
        res = pool.apply_async(build_bowtie2_index, (split_genome_abs_paths[i], ))
        result.append(res)
    pool.close()
    pool.join()
    bowtie2_index_abs_paths = []
    for res in result:
        bowtie2_index_abs_paths.append(res.get())

    # map the input reads file to split genomes
    pool = multiprocessing.Pool(processes = multiprocess_core)
    for i in xrange(len(bowtie2_index_abs_paths)):
        pool.apply_async(map_with_bowtie2, (bowtie2_index_abs_paths[i], real_data_path))
    pool.close()
    pool.join()

    # count reads in the bins of each genome
    sam_file_abs_paths = list_sam_file(genome_path)
    acc_micro_name_dict = acc_microbe_name(genome_path)

    pool = multiprocessing.Pool(processes = multiprocess_core)
    result = []
    for i in xrange(len(sam_file_abs_paths)):
        res = pool.apply_async(read_count, (sam_file_abs_paths[i], acc_micro_name_dict))
        result.append(res)
    pool.close()
    pool.join()
    genomes_distribution = []
    for res in result:
        genomes_distribution.append(res.get())    

    # save the genome distribution to a file named "distribution.txt"
    write_genome_distribution(genome_distribution_path, genomes_distribution)    
    
if __name__ == '__main__':
    main()
