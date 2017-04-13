"""
Function:     calculate the relationship between GC content and coverage bias 
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

def list_fastq_data(file_path):
    """
    Args:
        file_path:         str; path of the input fastq data, only ".fq" suffix is allowed                
    """
    file_path = file_path.rstrip('/') + '/'
    file_abs_paths = []
    file_names = os.listdir(file_path)
    for file_ in file_names:
        if file_.endswith(".fq"):
            file_name = file_path + file_
            file_abs_paths.append(file_name)
    return file_abs_paths

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

def fastq_genome(fastq_path, genome_path):
    """
    Args:
        fastq_path:     str; absolute paths of the input fastq data
        genome_path:    str; absolute paths of the input genomes
    Note:
        make sure the prefixes of the single genome sequencing fastq data and its corresponding reference 
        genome are the same, e.g., if the name of the fastq data is "single_genome.fq", its corresponding
        reference genome name should be "single_genome.fna", otherwise the program do not know the relationship
        between the single genome sequenicng data and its reference genome            
    """
    fastq_path_dict = {}
    for fq_path in fastq_path:
        fq_name_prefix = fq_path.split('/')[-1].rstrip(".fq")
        fastq_path_dict[fq_name_prefix] = fq_path
    fastq_genome_dict = {}
    for geo_path in genome_path:
        geo_name_prefix = geo_path.split('/')[-1].rstrip(".fna")
        if fastq_path_dict.has_key(geo_name_prefix):
            fastq_genome_dict[fastq_path_dict[geo_name_prefix]] = geo_path
    return fastq_genome_dict

def split_genome(fastq_path, fq_genome_dict, bin_size = 2000):
    """
    Args:
        fastq_path:        str; absolute paths of the input fastq data
        fq_genome_dict:    str; single genome sequencing data and its corresponding reference genome
        bin_size:        int; the size of each bin to be cut of the input genome    
    """

    def calculate_GC(sequence):
        """
        Args:
            sequence:        str; the sequence of a bin in the input genome
        """
        sequence = sequence.upper()
        seq_len = 0
        GC = 0
        for base in sequence:
            seq_len += 1
            if base == 'C' or base == 'G':
                GC += 1
        GC_content = np.int(np.round(GC / np.float(seq_len) * 100))
        return GC_content

    # read one genome
    annotation = ""
    genome = []
    genome_abs_path = fq_genome_dict[fastq_path]
    with open(genome_abs_path) as f_r:
        annotation = f_r.readline().strip()        # e.g., >U00096.3 Escherichia coli str. K-12 substr. MG1655, complete genome
        acc_num = annotation.split()[0].lstrip('>')
        for line in f_r:
            genome.extend(line.strip())
    genome = ''.join(genome)
    # split the genome using bin_size
    genome_len = len(genome)
    n_bins = np.int(np.ceil(genome_len / np.float(bin_size)))
    genome_path_prefix = '/'.join(genome_abs_path.split('/')[:-1]) + '/'
    genome_abs_path = genome_abs_path.split('/')[-1].rstrip(".fna") + "_split.fna"
    genome_abs_path = genome_path_prefix + genome_abs_path
    with open(genome_abs_path, 'w') as f_w:
        for i in xrange(n_bins - 1):
            seq = genome[i * bin_size:(i + 1) * bin_size]
            GC_content = calculate_GC(seq)
            annot = '>' + str(i) + '-' + "GC%s" % str(GC_content) + '-' + acc_num
            f_w.write(annot)
            f_w.write('\n')
            f_w.write(seq)
            f_w.write('\n')
        seq = genome[(n_bins - 1) * bin_size:]
        GC_content = calculate_GC(seq)
        annot = '>' + str(n_bins - 1) + '-' + "GC%s" % str(GC_content) + '-' + acc_num
        f_w.write(annot)
        f_w.write('\n')
        f_w.write(seq)
        f_w.write('\n')
    return fastq_path, genome_abs_path

def build_bowtie2_index(fastq_path, fq_genome_dict):
    """
    Args:
        fastq_path:        str; absolute paths of the input fastq data
        fq_genome_dict:    str; single genome sequencing data and its corresponding reference split genome
    """
    split_genome_abs_path = fq_genome_dict[fastq_path]
    split_genome_path = '/'.join(split_genome_abs_path.split('/')[:-1]) + '/'
    bowtie2_index = split_genome_abs_path.split('/')[-1].rstrip(".fna") + "_index"
    bowtie2_index = split_genome_path + bowtie2_index
    os.system("bowtie2-build %s %s" % (split_genome_abs_path, bowtie2_index))
    return fastq_path, bowtie2_index

def map_with_bowtie2(fastq_path, fq_index_dict):
    """
    Args:
        fastq_path:     str; absolute paths of the input fastq data            
        fq_index_dict:    str; single genome sequencing data and its corresponding built bowtie2 index
    """
    bowtie2_index_abs_path = fq_index_dict[fastq_path]
    bowtie2_index_path = '/'.join(bowtie2_index_abs_path.split('/')[:-1]) + '/'
    split_genome_name_prefix = bowtie2_index_abs_path.split('/')[-1].rstrip("_index")
    fastq_path_name_prefix = fastq_path.split('/')[-1].split('.')[0]
    sam_file_name = bowtie2_index_path + split_genome_name_prefix + '_' + fastq_path_name_prefix + ".sam"
    os.system("bowtie2 -x %s -U %s -q -S %s --no-unal" % (bowtie2_index_abs_path, fastq_path, sam_file_name))

def GC_read_count_relation(file_name):
    """
    Args:
        file_name:         str; the .sam file to be counted        
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
                    GC_content = np.int(bin_name.split('-')[1].lstrip("GC"))
                    read_count_dict[bin_num] = np.array([GC_content, 0])
                else:
                    continue
            else:
                bin_name = line.split()[2]
                bin_num = np.int(bin_name.split('-')[0])
                read_count_dict[bin_num][1] += 1
    GC_read_count_dict = {}
    for i in xrange(101):                             # 101 means 0%-100%, step is 1%
        GC_read_count_dict[i] = np.array([0, 0])
    for bin_num in read_count_dict:
        GC_content = read_count_dict[bin_num][0]
        read_count = read_count_dict[bin_num][1]
        GC_read_count_dict[GC_content] += np.array([1, read_count])
    for GC_content in GC_read_count_dict:
        GC_sample_number = GC_read_count_dict[GC_content][0]
        read_count = GC_read_count_dict[GC_content][1]
        if read_count == 0:
            GC_read_count_dict[GC_content] = 0
        else:
            GC_read_count_dict[GC_content] = read_count / np.float(GC_sample_number)
    GC_read_count_list = []
    for GC_content in GC_read_count_dict:
        GC_read_count_term = str(GC_content) + '-' + str(GC_read_count_dict[GC_content])
        GC_read_count_list.append(GC_read_count_term)
    GC_read_count_list = ':'.join(GC_read_count_list)
    return GC_read_count_list

def write_GC_read_count(file_path, GC_read_count):
    """
    Args:
        file_name:         str; the absolute folder path to save the file of relationship between GC content
                        and coverage bias
        GC_read_count:    str; a string denotes the relationship between GC content and coverage bias
    """
    file_path = file_path.rstrip('/') + '/'
    file_name = file_path + "GC_coverage.txt"
    with open(file_name, 'w') as f_w:
        f_w.write("# Relationship between GC content and coverage bias")
        f_w.write('\n')
        for term in GC_read_count:
            f_w.write(term)
            f_w.write('\n')
    print "GC coverage bias file saved in %s." % file_name

def main():
    usageReminder = ("usage: python cal_GC_coverage_relation_step1.py [--fq_path] fastq_data_folder\n"
                     "[--geo_path] ref_genome_folder [-o] output_path [-b] bin_size")
    parser = OptionParser(usage = usageReminder)
    parser.add_option("--fq_path", action = "store", type = "string", dest = "fastq_data_path", \
                      help = "Input folder path for all the single genome sequencing data sample.")
    parser.add_option("--geo_path", action = "store", type = "string", dest = "genome_path", \
                      help = ("Input folder path for all the reference genomes of each corresponding "
                                "single genome sequencing data sample."))
    parser.add_option("-o", action = "store", type = "string", dest = "output_path", \
                      help = "Output folder path for saving the GC-based coverage bias file.")
    parser.add_option("-b", action = "store", type = "int", dest = "bin_size", \
                      help = ("The sliding window size for calculating the relationship between GC content "
                                "and the mapped reads count. DeFault is 2000."))
    (opts, args) = parser.parse_args()
    if opts.bin_size == None:
        opts.bin_size = 2000

    fastq_data_path = opts.fastq_data_path
    genome_path = opts.genome_path
    GC_coverage_path = opts.output_path
    bin_size = opts.bin_size
    # fastq_data_path = "/data/ssliu/Simulation/SimExperimentPipe/simMicrobeEnv/genNonuniReads/single_genome_sequencing_data/temp_fastq"
    # genome_path = "/data/ssliu/Simulation/SimExperimentPipe/simMicrobeEnv/genNonuniReads/single_genome_sequencing_data/temp_genome"
    # GC_coverage_path = "/data/ssliu/Simulation/SimExperimentPipe/simMicrobeEnv/genNonuniReads/single_genome_sequencing_data/temp_fastq"
    # bin_size = 2000

    # build the link between the single genome sequencing data and its corresponding reference genome
    fastq_data_abs_paths = list_fastq_data(fastq_data_path)
    genome_abs_paths = list_genome(genome_path)
    fastq_genome_dict = fastq_genome(fastq_data_abs_paths, genome_abs_paths)
    
    multiprocess_core = multiprocessing.cpu_count() / 2    # defalut processing core number (half of the total cpu cores)

    # split genome
    pool = multiprocessing.Pool(processes = multiprocess_core)
    result = []
    for i in xrange(len(fastq_data_abs_paths)):
        res = pool.apply_async(split_genome, (fastq_data_abs_paths[i], fastq_genome_dict, bin_size))
        result.append(res)
    pool.close()
    pool.join()
    fastq_split_genome_dict = {}
    for res in result:
        fastq_path, genome_abs_path = res.get()
        fastq_split_genome_dict[fastq_path] = genome_abs_path

    # build bowtie2 index
    pool = multiprocessing.Pool(processes = multiprocess_core)
    result = []
    for i in xrange(len(fastq_data_abs_paths)):
        res = pool.apply_async(build_bowtie2_index, (fastq_data_abs_paths[i], fastq_split_genome_dict))
        result.append(res)
    pool.close()
    pool.join()
    fastq_bowtie2_index_dict = {}
    for res in result:
        fastq_path, bowtie2_index = res.get()
        fastq_bowtie2_index_dict[fastq_path] = bowtie2_index

    # map the input reads file to split genomes
    pool = multiprocessing.Pool(processes = multiprocess_core)
    for i in xrange(len(fastq_data_abs_paths)):
        pool.apply_async(map_with_bowtie2, (fastq_data_abs_paths[i], fastq_bowtie2_index_dict))
    pool.close()
    pool.join()

    # count reads in the bins of each genome
    sam_file_abs_paths = list_sam_file(genome_path)

    pool = multiprocessing.Pool(processes = multiprocess_core)
    result = []
    for i in xrange(len(sam_file_abs_paths)):
        res = pool.apply_async(GC_read_count_relation, (sam_file_abs_paths[i], ))
        result.append(res)
    pool.close()
    pool.join()
    GC_read_count = []
    for res in result:
        GC_read_count.append(res.get())

    # save the GC content to a file named "GC_coverage.txt"
    write_GC_read_count(GC_coverage_path, GC_read_count)

if __name__ == '__main__':
    main()
