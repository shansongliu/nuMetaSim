"""
Function:     generate strains of the input genome  
"""

import os
import re
import numpy as np
import multiprocessing
from optparse import OptionParser

def argvParse():
    usageReminder = ("usage: python gen_genome_strain.py [--geo_va] genome_variation_file [--geo] genome_file\n"
                     "[-o] output_folder [-m] mode_flag [-t] variation_rates [-n] genome_number\n"
                     "[--args] normal_distri_args [-v] var_seed_flag [-r] random_seed")
    parser = OptionParser(usage = usageReminder)
    parser.add_option("--geo_va", action = "store", type = "string", dest = "genome_variation_file", \
                      help = "Genome variation file generated by \"cal_genome_variation.py.\"")
    parser.add_option("--geo", action = "store", type = "string", dest = "genome_file", \
                      help = "Genome file used for \"relative strain\" generation.")
    parser.add_option("-o", action = "store", type = "string", dest = "output_folder", \
                        help = "Output folder for saving generated \"relative strain\" genomes.")
    parser.add_option("-m", action = "store", type = "int", dest = "run_mode_flag", \
                      help = ("A 0-1 flag to determine whether the user inputs all the variation rates "
                                "or only inputs the genome number to be generated and the arguments of "
                                "normal distribution (mean, variation). 0 means the user provides a list of "
                                "variation rates. The input format is like \"0.1,0.2,0.3,0.4,0.5\", then the "
                                "program will use the list of rates to generate a series of varied genomes. "
                                "The genome number to be generated in this mode is determined by the length of "
                                "the input variation rates. 1 means the program will use normal distribution to "
                                "generate variation rates. But the genome number to be generated and the arguments "
                                "(mean and variation, input format is \"0,1\") of the normal distribution should "
                                "be input by the user. Note that the maximum variation rate cannot be higher than "
                                "0.5, the program will check all the variation rates to ensure their validity. If "
                                "some of the rates are higher than 0.5, they will be forced to 0.5."))
    parser.add_option("-t", action = "store", type = "string", dest = "variation_rates", \
                      help = "Variation rates input by the user, please refer the input format in \"-m\" parameter.")
    parser.add_option("-n", action = "store", type = "int", dest = "genome_num", \
                      help = "Genome number to be generated, only works when the \"-m\" parameter is set to be 1.")
    parser.add_option("--args", action = "store", type = "string", dest = "normal_distri_args", \
                      help = ("Normal distribution arguments, including mean and variation. only works when the \"-m\" "
                                "parameter is set to be 1. Please refer the input format in the \"-m\" parameter. "
                                "Default is \"0.1,0.1\"."))
    parser.add_option("-v", action = "store", type = "int", dest = "var_seed_flag", \
                      help = ("A 0-1 flag to control whether to use a variable random seed each running time. 0 means "
                              "fixed random seed (should be set by the user), 1 means variable random seed. Default is 0."))
    parser.add_option("-r", action = "store", type = "int", dest = "random_seed", \
                      help = "Randon seed, only works when the \"var_seed_flag\" is set be 0. Default random_seed is 0.")
    (opts, args) = parser.parse_args()
    return opts    

def argvCheck(opts):
    # check run mode flag
    if opts.run_mode_flag == None or opts.run_mode_flag == 0:
        opts.run_mode_flag = 0
    # check variation rates
    if opts.run_mode_flag == 0 and opts.variation_rates == None:
        print "The running mode is 0, please provide the variation rates. Program will exit."
        os._exit(1)
    elif opts.run_mode_flag == 0 and opts.variation_rates != None:
        p = re.compile(r'-?[0-9]+\.?[0-9]*')
        tmp_list = p.findall(opts.variation_rates)
        opts.variation_rates = np.array([float(term) for term in tmp_list])
    elif opts.run_mode_flag == 1 and opts.genome_num == None:
        print "The running mode is 1, please provide the genome number. Program will exit."
        os._exit(1)
    elif opts.run_mode_flag == 1 and opts.genome_num != None:
        if opts.normal_distri_args == None:
            opts.normal_distri_args = np.array([0.1, 0.1])
        elif opts.normal_distri_args != None:
            p = re.compile(r'-?[0-9]+\.?[0-9]*')
            tmp_list = p.findall(opts.normal_distri_args)
            opts.normal_distri_args = np.array([float(term) for term in tmp_list])
    # check variable seed flag
    if opts.var_seed_flag == None or opts.var_seed_flag == 0:
        opts.var_seed_flag = 0
        if opts.random_seed == None:
            opts.random_seed = 0
    return opts

def generate_variation_rate(genome_num, mean, std):
    """
    Args:
        genome_num:        int; the set number of genomes to be generated
        mean:            float; the set mean for generating variation rate
        std:            float; the set standard deviation for generating variation rate
    """
    variation_rate_threshold = 0.5
    rates = np.random.normal(mean, std, genome_num)
    rates = np.abs(rates)
    for i in xrange(len(rates)):
        if rates[i] > variation_rate_threshold:
            rates[i] = variation_rate_threshold
    return rates

def read_genome_variation(file_name):
    """
    Args:
        file_name:        str; genome variation file "genome_variation.txt" generated by "cal_genome_variation.py"
                        or a user-defined genome variation file in accordance with the format of "genome_variation.txt"
    """
    with open(file_name) as f_r:
        f_r.readline()
        # relative proportion of substitution/deletion/insertion
        sub_del_ins_propor = f_r.readline().strip()
        sub_del_ins_propor = sub_del_ins_propor.split('=')[-1].split(':')
        sub_del_ins_propor = np.array([np.float(propor) for propor in sub_del_ins_propor]) # sub:del:ins
        # base substitution probability
        sub_base_propor = f_r.readline().strip()
        sub_base_propor = sub_base_propor.split('=')[-1].split(':')
        sub_base_propor = np.array([np.float(propor) for propor in sub_base_propor]) # AA:AG:AC:AT:GA:GG:GC:GT:CA:CG:CC:CT:TA:TG:TC:TT
        sub_base_propor = sub_base_propor.reshape((4, 4)) #   A G C T
                                                          # A
                                                          # G
                                                          # C
                                                          # T
        # substitution length distribution
        sub_length_distribution = f_r.readline().strip()
        sub_length_distribution = sub_length_distribution.split('=')[-1].split(':')
        sub_length_distribution = np.array([np.float(term.split('_')[-1]) for term in sub_length_distribution])
        # deletion length distribution
        del_length_distribution = f_r.readline().strip()
        del_length_distribution = del_length_distribution.split('=')[-1].split(':')
        del_length_distribution = np.array([np.float(term.split('_')[-1]) for term in del_length_distribution])
        # base insertion probability
        ins_base_propor = f_r.readline().strip()
        ins_base_propor = ins_base_propor.split('=')[-1].split(':')
        ins_base_propor = np.array([np.float(propor) for propor in ins_base_propor]) # A:G:C:T
        # insertion length distribution
        ins_length_distribution = f_r.readline().strip()
        ins_length_distribution = ins_length_distribution.split('=')[-1].split(':')
        ins_length_distribution = np.array([np.float(term.split('_')[-1]) for term in ins_length_distribution])
    return sub_del_ins_propor, sub_base_propor, sub_length_distribution, del_length_distribution, ins_base_propor, ins_length_distribution 

def read_genome(file_name):
    """
    Args:
        file_name:      str; the genome file, it should only contains a single complete sequence
                        the file format is .fna (using .fa for short), which contains:
                            >U00096.3 Escherichia coli str. K-12 substr. MG1655, complete genome
                            ATCGAGTCAGTCAGTCGATCG...
                            CAGTCACGATAGTATCGGTCG...
                            GTCAGCAAGTCTCGATCGAGT...                
    """
    annotation = ""
    genome = []
    with open(file_name) as f_r:
        annotation = f_r.readline().strip() # e.g., >U00096.3 Escherichia coli str. K-12 substr. MG1655, complete genome
        for line in f_r:
            genome.extend(line.strip())
    genome = ''.join(genome)
    return annotation, genome

def modify_genome(genome, annotation, sdi_propor, sub_base_propor, sub_len_distri,
                  del_len_distri, ins_base_propor, ins_len_distri, count, file_name,
                  file_path_write, variation_rate = 0.05):
    """
    Args:
        genome:            str; a string denotes the genome, like 'ATCGACTGAA...'
        annotation:        str; an annotation string extracted from the genome file
        sdi_propor:        np.array; relative proportion of substitution/deletion/insertion
        sub_base_propor:np.array; base substitution probability
        sub_len_distri:    np.array; substitution length distribution
        del_len_distri:    np.array; deletion length distribution        
        ins_base_propor:np.array; base insertion probability
        ins_len_distri:    np.array; insertion length distribution
        count:            int; the current serial number of the modified genome
        file_name:        str; the input genome file's absolute path
        file_path_write:str; the output directory for modified genomes
        variation_rate: float; genome variation rate, say 0.1, then about 10 percent bases of the input genome 
                        will be modified (subtituted/inserted/deleted), default is 0.05, no larger than 0.3
    """
    def substitution(sub_base_propor, original):
        if original == 'A':
            ori = 0
        elif original == 'G':
            ori = 1
        elif original == 'C':
            ori = 2
        else:
            ori = 3 # base 'T'
        dice = np.random.uniform(0, 1)
        if dice < sub_base_propor[ori][0]:
            return 'A' # sub to 'A'
        elif dice < sub_base_propor[ori][0] + sub_base_propor[ori][1]:
            return 'G' # sub to 'G'
        elif dice < sub_base_propor[ori][0] + sub_base_propor[ori][1] + sub_base_propor[ori][2]:
            return 'C' # sub to 'C'
        else:
            return 'T' # sub to 'T'

    def insertion(ins_base_propor):
        dice = np.random.uniform(0, 1)
        if dice < ins_base_propor[0]:
            return 'A' # sub to 'A'
        elif dice < ins_base_propor[0] + ins_base_propor[1]:
            return 'G' # sub to 'G'
        elif dice < ins_base_propor[0] + ins_base_propor[1] + ins_base_propor[2]:
            return 'C' # sub to 'C'
        else:
            return 'T' # sub to 'T'    

    try:
        import intervaltree
    except ImportError:
        os.system("pip install intervaltree", shell = True)
    import intervaltree

    SDI = ['S', 'D', 'I']
    variation_rate_threshold = 0.5
    if variation_rate > variation_rate_threshold:
        variation_rate = variation_rate_threshold
        print "The set variation rate is larger than the allowed maximum (0.5), now the program will use 0.5 instead"
    genome_len = len(genome)
    genome = list(genome)

    # calculate the modified base number of each variation
    tot_modified_bases = genome_len * variation_rate
    sub_del_ins_bases = sdi_propor * tot_modified_bases
    sub_del_ins_distri = [sub_len_distri, del_len_distri, ins_len_distri]
    sub_del_ins_lengths = []
    sub_del_ins_len_avg = []
    length_max = []
    for i in xrange(len(sub_del_ins_distri)):
        lengths = np.arange(len(sub_del_ins_distri[i])) + 1
        length_max.append(lengths[-1])
        sub_del_ins_lengths.append(lengths)
        len_avg = np.sum(lengths * sub_del_ins_distri[i])
        sub_del_ins_len_avg.append(len_avg)
    length_max = np.max(length_max)
    sub_del_ins_len_avg = np.array(sub_del_ins_len_avg)
    sub_del_ins_sample_points_num = sub_del_ins_bases / sub_del_ins_len_avg
    sub_del_ins_sample_points_num = map(lambda x: np.int(np.round(x)), sub_del_ins_sample_points_num)    
    # generate sample points needed to be modified on the input genome
    sample_num = np.sum(sub_del_ins_sample_points_num)
    sample_points = np.random.uniform(0, 1, sample_num) * (genome_len - length_max)
    sample_points = np.array(sample_points, dtype = np.int)
    # label the positions to be modified with modification type and length
    labels = []
    for i in xrange(len(sub_del_ins_sample_points_num)):
        labels.extend(list(SDI[i] * sub_del_ins_sample_points_num[i]))
    # generate lengths using substitution/deletion/insertion length distribution
    length_samples = []
    interval_tree_dict = {}
    for i in xrange(len(sub_del_ins_distri)):
        len_distri_cumsum = np.cumsum(sub_del_ins_distri[i])
        interval_tree = intervaltree.IntervalTree()
        interval_tree.addi(0, sub_del_ins_distri[i][0], sub_del_ins_lengths[i][0])
        for j in xrange(1, len(len_distri_cumsum)):
            interval_tree.addi(len_distri_cumsum[j - 1], len_distri_cumsum[j], sub_del_ins_lengths[i][j])
        interval_tree_dict[SDI[i]] = interval_tree
        length_samples_raw = np.random.uniform(0, 1, sub_del_ins_sample_points_num[i])
        for sample in length_samples_raw:
            len_sample = list(interval_tree.search(sample))[0].data
            length_samples.append(len_sample)
    geo_vary_points_label = zip(labels, length_samples)
    geo_vary_points = zip(sample_points, geo_vary_points_label)
    geo_vary_points = sorted(geo_vary_points, key = lambda x: x[0])
    # modify genome
    for i in xrange(len(geo_vary_points) - 1, -1, -1):
        geo_ind, vary_type, vary_len = geo_vary_points[i][0], geo_vary_points[i][1][0], geo_vary_points[i][1][1]
        if vary_type == 'S':
            for j in xrange(geo_ind, geo_ind + vary_len):
                genome[j] = substitution(sub_base_propor, genome[geo_ind])
        elif vary_type == 'D':
            for j in xrange(geo_ind, geo_ind + vary_len):
                genome[j] = 'D'
        elif vary_type == 'I':
            for j in xrange(geo_ind, geo_ind + vary_len):
                inserted_base = insertion(ins_base_propor)
                genome.insert(j, inserted_base)
    genome = ''.join(genome)
    genome = re.sub('D', '', genome)
    # write genome to file
    write_genome(genome, annotation, count, file_name, file_path_write)

def write_genome(genome, annotation, count, file_name, file_path_write):
    """
    Args:
        genome:            str; modified genome given the input genome and variation rate
        annotation:        str; an annotation string extracted from the genome file
        count:            int; the current serial number of the modified genome
        file_name:        str; the input genome file's absolute path
        file_path_write:str; the output directory for modified genomes
    """
    file_name = file_name.split('/')[-1].rstrip(".fna")
    file_path_write = file_path_write.rstrip('/') + '/'
    file_name = file_path_write + file_name + "_modified_" + str(count) + ".fna"
    if ',' in annotation:
        annotation = annotation.split(',')[0]
    annotation += " relative strain " + str(count)

    genome_len = len(genome)
    line_len = 80
    line_num = genome_len / line_len
    with open(file_name, 'w') as f_w:
        f_w.write(annotation)
        f_w.write('\n')
        for i in xrange(line_num):
            write_line = genome[i * line_len:(i + 1) * line_len]
            f_w.write(write_line)
            f_w.write('\n')
        if genome_len - line_num * line_len > 0:
            write_line = genome[line_num * line_len:]
            f_w.write(write_line)
            f_w.write('\n')

def main():
    opts = argvParse()
    opts = argvCheck(opts)

    genome_variation_file = opts.genome_variation_file
    genome_file = opts.genome_file
    file_path_write = opts.output_folder
    mode_flag = opts.run_mode_flag
    variation_rates = opts.variation_rates
    genome_num = opts.genome_num
    normal_distri_args = opts.normal_distri_args
    mean, std = 0.0, 0.0
    if mode_flag == 1:
        mean, std = normal_distri_args[0], normal_distri_args[1]
    var_seed_flag = opts.var_seed_flag
    random_seed = opts.random_seed

    if mode_flag == 1:
        variation_rates = generate_variation_rate(genome_num, mean, std)
    # set a global fixed random seed if needed
    if var_seed_flag == 0:
        np.random.seed(random_seed)
    else:
        import time
        np.random.seed(np.int(time.time()))    

    sdi_propor, sub_base_propor, sub_len_distri, del_len_distri, ins_base_propor, ins_len_distri = read_genome_variation(genome_variation_file)
    annotation, genome = read_genome(genome_file)
    # generate modified genomes
    multiprocess_core = multiprocessing.cpu_count() / 2    # defalut processing core number (half of the total cpu cores)
    pool = multiprocessing.Pool(processes = multiprocess_core)
    count = 0
    for variation_rate in variation_rates:
        pool.apply_async(modify_genome, (genome, annotation, sdi_propor, sub_base_propor, sub_len_distri, 
                                         del_len_distri, ins_base_propor, ins_len_distri, count, genome_file,
                                           file_path_write, variation_rate))
        count += 1
    pool.close()
    pool.join()

if __name__ == '__main__':
    main()
