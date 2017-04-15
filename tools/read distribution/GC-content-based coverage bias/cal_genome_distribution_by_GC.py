"""
Function:     calculate the genome distribution by GC content 
"""

import os
import numpy as np
import multiprocessing
from optparse import OptionParser

def read_polyfit_coefficient(file_name):
    """
    Args:
        file_name:         str; absolute path of "polyfit_coefficient.txt"                
    """
    with open(file_name) as f_r:
        f_r.readline()
        coefs = f_r.readline().strip().split()
        coefs = [np.float(coef) for coef in coefs]
        coefs = np.array(coefs)
    return coefs

def read_genome(file_name):
    """
    Args:
        file_name:      str; the genome file, it should only contains a single complete sequence              
    """
    annotation = ""
    genome = []
    with open(file_name) as f_r:
        annotation = f_r.readline().strip()        # e.g., >U00096.3 Escherichia coli str. K-12 substr. MG1655, complete genome
        for line in f_r:
            genome.extend(line.strip())
    genome = ''.join(genome)
    return annotation, genome

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

def cal_genome_distri_by_GC(file_name, polyfit_coef, bin_size = 2000):
    """
    Args:
        file_name:         str; the genome file, it should only contains a single complete sequence
        polyfit_coef:    np.array; polynomial fitting coefficient
        bin_size:       int; the size of the bins in the input genome, in order to seperate the genome
                          into several bins and assign different probabilities of each bin, to imitate a
                          non-uniform probability distribution to extract the reads from the input genome        
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

    annotation, genome = read_genome(file_name)
    microbe_name = annotation.lstrip('>').split(' ', 1)[-1]
    if ',' in microbe_name:
        microbe_name = microbe_name.split(',')[0]
    genome_len = len(genome)
    n_bins = np.int(np.ceil(genome_len / np.float(bin_size)))
    GC_contents = np.zeros(n_bins)
    for i in xrange(n_bins - 1):
        seq = genome[i * bin_size:(i + 1) * bin_size]
        GC_content = calculate_GC(seq)
        GC_contents[i] = GC_content
    seq = genome[(n_bins - 1) * bin_size:]
    GC_content = calculate_GC(seq)
    GC_contents[n_bins - 1] = GC_content
    
    # transform GC to proportion
    poly1d_coef = np.poly1d(polyfit_coef)
    fitted_propor = poly1d_coef(GC_contents)
    fitted_propor *= 100000        # scale the proportion to larger value
    pseudo_count = 0.01            # to avoid 0 count genome bin, we add 0.01 to this kind of bin
    fitted_propor[fitted_propor <= 0] = pseudo_count
    fitted_propor_list = []
    for num in fitted_propor:
        fitted_propor_list.append(str(num))
    fitted_propor_list = ':'.join(fitted_propor_list)

    return (microbe_name, fitted_propor_list)

def write_genome_distribution(file_name, geo_distri):
    """
    Args:
        file_name:         str; absolute path of "polyfit_coefficient.txt"    
        geo_distri:        list; a list contains each calculated genome's distribution by counting reads 
                        which were mapped to bins in the genome            
    """
    file_path = ""
    if '/' in file_name:
        file_path = '/'.join(file_name.split('/')[:-1]) + '/'
    else:
        file_path = "./"
    file_name = file_path + "distribution.txt"
    with open(file_name, 'w') as f_w:
        for term in geo_distri:
            microbe_name, distribution = term[0], term[1]
            write_line = microbe_name + '\t' + distribution
            f_w.write(write_line)
            f_w.write('\n')
    print "Genome distribution file calculated by GC content saved in %s." % file_name

def main():
    usageReminder = ("usage: python cal_genome_distribution_by_GC.py [-f] GC_fitting_coef [--geo_path] ref_genome_folder\n"
                     "[-b] bin_size")
    parser = OptionParser(usage = usageReminder)
    parser.add_option("-f", action = "store", type = "string", dest = "GC_fitting_coef", \
                      help = "The fitting coefficient file calculated by \"polyfit_GC_coverage.py.\"")
    parser.add_option("--geo_path", action = "store", type = "string", dest = "ref_genome_folder", \
                      help = "The folder path of all the reference genomes for GC-based genome distribution calculation.")
    parser.add_option("-b", action = "store", type = "int", dest = "bin_size", \
                      help = ("The sliding window size for calculating the reads proportion to be generated in a local "
                              "genomic region of size \"bin_size\". DeFault is 2000."))
    (opts, args) = parser.parse_args()
    if opts.bin_size == None:
        opts.bin_size = 2000
    
    GC_coverage_path = opts.GC_fitting_coef
    genome_path = opts.ref_genome_folder
    bin_size = opts.bin_size

    # GC_coverage_path = "/data/ssliu/Simulation/SimExperimentPipe/simMicrobeEnv/genNonuniReads/real_seq_data/polyfit_coefficient.txt"
    # genome_path = "/data/ssliu/Simulation/SimExperimentPipe/simMicrobeEnv/genNonuniReads/microbes"
    # bin_size = 2000

    multiprocess_core = multiprocessing.cpu_count() / 2    # defalut processing core number (half of the total cpu cores)

    # calculate distribution of each genome by GC content
    coef = read_polyfit_coefficient(GC_coverage_path)
    genome_abs_paths = list_genome(genome_path)
    pool = multiprocessing.Pool(processes = multiprocess_core)
    result = [] 
    for i in xrange(len(genome_abs_paths)):
        res = pool.apply_async(cal_genome_distri_by_GC, (genome_abs_paths[i], coef, bin_size))
        result.append(res)
    pool.close()
    pool.join()
    genomes_distribution = []
    for res in result:
        genomes_distribution.append(res.get())

    # save the calculated distribution file to "distribution.txt"
    write_genome_distribution(GC_coverage_path, genomes_distribution)

if __name__ == '__main__':
    main()
