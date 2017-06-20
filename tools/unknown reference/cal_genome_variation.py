"""
Function:     calculate the insertion/deletion/substitution  
"""

import os
import numpy as np
from optparse import OptionParser

def align_ref_qry(genome_ref, genome_qry):
    cur_work_dir = os.path.abspath(os.curdir)
    dir_to_change = ""
    if '/' in genome_ref:
        dir_to_change = '/'.join(genome_ref.split('/')[:-1])
    else:
        dir_to_change = "./"
    os.chdir(dir_to_change)
    # build reference genome index
    try:
        os.system("makeblastdb -in %s -dbtype nucl" % genome_ref)
    except:
        print "Please make sure the required program \"makeblastdb\" is installed in the system path. Program will exit."
        os._exit(1)
    # align query genome to reference genome
    dir_to_change = ""
    if '/' in genome_qry:
        dir_to_change = '/'.join(genome_qry.split('/')[:-1])
    else:
        dir_to_change = "./"
    os.chdir(dir_to_change)
    dir_to_change = dir_to_change.rstrip('/') + '/'    
    alignment_file = dir_to_change + "alignment.out"
    outfmt = "6 qseqid sseqid pident length mismatch gapopen gaps qstart qend sstart send qseq sseq evalue bitscore"
    try:
        os.system("blastn -db %s -query %s -out %s -outfmt \"%s\" -evalue 1e-5" % (genome_ref, genome_qry, alignment_file, outfmt))
        print ("Alignment file bwtween two genomes saved in %s.") % alignment_file        
    except:
        print "Please make sure the required program \"blastn\" is installed in the system path. Program will exit."
        os._exit(1)
    os.chdir(cur_work_dir)    
    return alignment_file

def read_genome_align_res(file_name):
    """
    Args:
        file_name:      str; the absolute path of the alignment file between a query genome
                        and a reference genome                
    """
    querys, refs = [], []
    with open(file_name) as f_r:
        for line in f_r:
            line = line.strip()
            terms = line.split()
            query, ref = terms[11], terms[12]
            querys.append(query)
            refs.append(ref)
    return querys, refs

def cal_variation(query, ref, sub_base_dict, sub_length_dict, del_length_dict, ins_base_dict, ins_length_dict):
    """
    Args:
        query:          str; aligned query sequence
        ref:            str; aligned reference sequence
        sub_base_dict:  dict; base substitution distribution        
        sub_length_dict:dict; substitution length distribution
        del_length_dict:dict; deletion length distribution
        ins_base_dict:  dict; base insertion distribution
        ins_length_dict:dict; insertion length distribution
    """
    BASE = "AGCT"
    qry_ref_len = len(query)
    qry_ind, ref_ind = 0, 0
    while qry_ind < qry_ref_len and ref_ind < qry_ref_len:
        if query[qry_ind] == ref[ref_ind]:            # no operation on the matched base
            qry_ind += 1
            ref_ind += 1
            continue
        else:                                        # insertion or deletion or substitution
            if query[qry_ind] in BASE:                # insertion or substitution
                if ref[ref_ind] in BASE:            # substitution
                    sub_orientation = ref[ref_ind] + "->" + query[qry_ind]
                    sub_base_dict[sub_orientation] += 1
                    qry_tmp_ind = qry_ind + 1
                    ref_tmp_ind = ref_ind + 1
                    sub_length = 1
                    while qry_tmp_ind < qry_ref_len and ref_tmp_ind < qry_ref_len:
                        if (query[qry_tmp_ind] != ref[ref_tmp_ind]) and (query[qry_tmp_ind] in BASE) and (ref[ref_tmp_ind] in BASE):
                            sub_orientation = ref[ref_tmp_ind] + "->" + query[qry_tmp_ind]
                            sub_base_dict[sub_orientation] += 1
                            sub_length += 1
                        else:
                            break
                        qry_tmp_ind += 1
                        ref_tmp_ind += 1
                    if not sub_length_dict.has_key(sub_length):
                        sub_length_dict[sub_length] = 1
                    else:
                        sub_length_dict[sub_length] += 1
                    qry_ind = qry_tmp_ind
                    ref_ind = ref_tmp_ind    
                else:                                # insertion
                    ins_base_dict[query[qry_ind]] += 1
                    qry_tmp_ind = qry_ind + 1
                    ref_tmp_ind = ref_ind + 1
                    ins_length = 1
                    while qry_tmp_ind < qry_ref_len and ref_tmp_ind < qry_ref_len:
                        if (query[qry_tmp_ind] != ref[ref_tmp_ind]) and (query[qry_tmp_ind] in BASE) and (ref[ref_tmp_ind] == '-'):
                            ins_base_dict[query[qry_tmp_ind]] += 1
                            ins_length += 1
                        else:
                            break        
                        qry_tmp_ind += 1
                        ref_tmp_ind += 1
                    if not ins_length_dict.has_key(ins_length):
                        ins_length_dict[ins_length] = 1
                    else:
                        ins_length_dict[ins_length] += 1
                    qry_ind = qry_tmp_ind
                    ref_ind = ref_tmp_ind
            else:                                    # deletion
                qry_tmp_ind = qry_ind + 1
                ref_tmp_ind = ref_ind + 1
                del_length = 1
                while qry_tmp_ind < qry_ref_len and ref_tmp_ind < qry_ref_len:
                    if (query[qry_tmp_ind] != ref[ref_tmp_ind]) and (query[qry_tmp_ind] == '-') and (ref[ref_tmp_ind] in BASE):
                        del_length += 1
                    else:
                        break        
                    qry_tmp_ind += 1
                    ref_tmp_ind += 1
                if not del_length_dict.has_key(del_length):
                    del_length_dict[del_length] = 1
                else:
                    del_length_dict[del_length] += 1
                qry_ind = qry_tmp_ind
                ref_ind = ref_tmp_ind

def write_genome_variation(sub_base_dict, sub_length_dict, del_length_dict, ins_base_dict, ins_length_dict, file_path):
    """
    Args:
        sub_base_dict:  dict; base substitution distribution        
        sub_length_dict:dict; substitution length distribution
        del_length_dict:dict; deletion length distribution
        ins_base_dict:  dict; base insertion distribution
        ins_length_dict:dict; insertion length distribution
        file_path:      str; the absolute path write the genome variation
    """
    # smoothing (use very small value to act as a pseudo count)
    pseudo_count = 1e-6
    # substitution smoothing
    for key in sub_base_dict.keys():
        if sub_base_dict[key] == 0:
            sub_base_dict[key] = 1e-6
    sub_lengths = sub_length_dict.keys()
    sub_length_min, sub_length_max = np.min(sub_lengths), np.max(sub_lengths)
    for i in xrange(sub_length_min, sub_length_max + 1):
        if (not sub_length_dict.has_key(i)) or (sub_length_dict[i] == 0):
            sub_length_dict[i] = pseudo_count
    # deletion smoothing
    del_lengths = del_length_dict.keys()
    del_length_min, del_length_max = np.min(del_lengths), np.max(del_lengths)
    for i in xrange(del_length_min, del_length_max + 1):
        if (not del_length_dict.has_key(i)) or (del_length_dict[i] == 0):
            del_length_dict[i] = pseudo_count
    # insertion smoothing
    for key in ins_base_dict.keys():
        if ins_base_dict[key] == 0:
            ins_base_dict[key] = 1e-6
    ins_lengths = ins_length_dict.keys()
    ins_length_min, ins_length_max = np.min(ins_lengths), np.max(ins_lengths)
    for i in xrange(ins_length_min, ins_length_max + 1):
        if (not ins_length_dict.has_key(i)) or (ins_length_dict[i] == 0):
            ins_length_dict[i] = pseudo_count

    # normalization
    # relative proportion of substitution/deletion/insertion
    sub_points = np.sum(sub_length_dict.values())
    del_points = np.sum(del_length_dict.values())
    ins_points = np.sum(ins_length_dict.values())
    sub_del_ins_tot_points = sub_points + del_points + ins_points
    sub_propor = sub_points / np.float(sub_del_ins_tot_points)
    del_propor = del_points / np.float(sub_del_ins_tot_points)
    ins_propor = ins_points / np.float(sub_del_ins_tot_points)
    sub_del_ins_relative_propor = "sub:del:ins=%f:%f:%f" % (sub_propor, del_propor, ins_propor)
    # base substitution probability 
    sub_tot_AtoX = sub_base_dict["A->G"] + sub_base_dict["A->C"] + sub_base_dict["A->T"]
    sub_tot_GtoX = sub_base_dict["G->A"] + sub_base_dict["G->C"] + sub_base_dict["G->T"]
    sub_tot_CtoX = sub_base_dict["C->A"] + sub_base_dict["C->G"] + sub_base_dict["C->T"]
    sub_tot_TtoX = sub_base_dict["T->A"] + sub_base_dict["T->G"] + sub_base_dict["T->C"]
    sub_base_dict["A->G"] /= np.float(sub_tot_AtoX)
    sub_base_dict["A->C"] /= np.float(sub_tot_AtoX)
    sub_base_dict["A->T"] /= np.float(sub_tot_AtoX)
    sub_base_dict["G->A"] /= np.float(sub_tot_GtoX)
    sub_base_dict["G->C"] /= np.float(sub_tot_GtoX)
    sub_base_dict["G->T"] /= np.float(sub_tot_GtoX)
    sub_base_dict["C->A"] /= np.float(sub_tot_CtoX)
    sub_base_dict["C->G"] /= np.float(sub_tot_CtoX)
    sub_base_dict["C->T"] /= np.float(sub_tot_CtoX)
    sub_base_dict["T->A"] /= np.float(sub_tot_TtoX)
    sub_base_dict["T->G"] /= np.float(sub_tot_TtoX)
    sub_base_dict["T->C"] /= np.float(sub_tot_TtoX)
    sub_base_relative_propor = "AA:AG:AC:AT:GA:GG:GC:GT:CA:CG:CC:CT:TA:TG:TC:TT=%f:%f:%f:%f:%f:%f:%f:%f:%f:%f:%f:%f:%f:%f:%f:%f" % \
                               (0.0, sub_base_dict["A->G"], sub_base_dict["A->C"], sub_base_dict["A->T"],
                                   sub_base_dict["G->A"], 0.0, sub_base_dict["G->C"], sub_base_dict["G->T"],
                                   sub_base_dict["C->A"], sub_base_dict["C->G"], 0.0, sub_base_dict["C->T"],
                                   sub_base_dict["T->A"], sub_base_dict["T->G"], sub_base_dict["T->C"], 0.0)
    # substitution length distribution normalization
    sub_length_normalizer = np.sum(sub_length_dict.values())
    for key in sub_length_dict.keys():
        sub_length_dict[key] /= np.float(sub_length_normalizer)
    sub_length_distribution = []
    for key in sub_length_dict.keys():
        sub_length_distribution.append(str(key) + '_' + str(sub_length_dict[key]))
    sub_length_distribution = ':'.join(sub_length_distribution)
    sub_length_distribution = "sub_length_distribution=" + sub_length_distribution
    # deletion length distribution normalization
    del_length_normalizer = np.sum(del_length_dict.values())
    for key in del_length_dict.keys():
        del_length_dict[key] /= np.float(del_length_normalizer)
    del_length_distribution = []
    for key in del_length_dict.keys():
        del_length_distribution.append(str(key) + '_' + str(del_length_dict[key]))
    del_length_distribution = ':'.join(del_length_distribution)
    del_length_distribution = "del_length_distribution=" + del_length_distribution
    # base insertion probability
    ins_base_normalizer = np.sum(ins_base_dict.values())
    for key in ins_base_dict.keys():
        ins_base_dict[key] /= np.float(ins_base_normalizer)
    ins_base_distribution = [str(value) for value in ins_base_dict.values()]
    ins_base_distribution = ':'.join(ins_base_distribution)
    ins_base_distribution = "ins_base_distribution=" + ins_base_distribution
    # insertion length distribution normalization
    ins_length_normalizer = np.sum(ins_length_dict.values())
    for key in ins_length_dict.keys():
        ins_length_dict[key] /= np.float(ins_length_normalizer)
    ins_length_distribution = []
    for key in ins_length_dict.keys():
        ins_length_distribution.append(str(key) + '_' + str(ins_length_dict[key]))
    ins_length_distribution = ':'.join(ins_length_distribution)
    ins_length_distribution = "ins_length_distribution=" + ins_length_distribution

    # write to file
    file_path = file_path.rstrip('/') + '/'
    file_name = file_path + "genome_variation.txt"
    with open(file_name, 'w') as f_w:
        f_w.write("#genome variation information")
        f_w.write('\n')
        f_w.write(sub_del_ins_relative_propor)
        f_w.write('\n')
        f_w.write(sub_base_relative_propor)
        f_w.write('\n')
        f_w.write(sub_length_distribution)
        f_w.write('\n')
        f_w.write(del_length_distribution)
        f_w.write('\n')
        f_w.write(ins_base_distribution)
        f_w.write('\n')
        f_w.write(ins_length_distribution)
        f_w.write('\n')
    print ("Genome variation file (prepared for \"gen_genome_strain.py\") saved in %s.") % file_name    

def main():
    usageReminder = ("usage: python cal_genome_variation.py [--ref] ref_genome [--qry] qry_genome [-d] file_save_path")
    parser = OptionParser(usage = usageReminder)
    parser.add_option("--ref", action = "store", type = "string", dest = "ref_genome_path", \
                      help = "The reference genome for genome variation calculation.")
    parser.add_option("--qry", action = "store", type = "string", dest = "qry_genome_path", \
                      help = "The query genome for genome variation calculation.")
    parser.add_option("-d", action = "store", type = "string", dest = "genome_variation_folder", \
                      help = ("Path to save the genome variation file between two genomes. The genome variation file "
                              "contains six parts:\n1.relative proportion of substitution/deletion/insertion\n"
                              "2.relative proportion of different substitution possibilities\n"
                              "3.substitution length distribution\n"
                              "4.deletion length distribution\n"
                              "5.relative proportion of different base insertion possibilities\n"
                              "6.insertion length distribution"))
    (opts, args) = parser.parse_args()

    sub_base_dict = {
        "A->G": 0, "A->C": 0, "A->T": 0, 
        "G->A": 0, "G->C": 0, "G->T": 0,
        "C->A": 0, "C->G": 0, "C->T": 0,
        "T->A": 0, "T->G": 0, "T->C": 0
    }
    sub_length_dict = {}
    del_length_dict = {}
    ins_base_dict = {"A": 0, "G": 0, "C": 0, "T":0}
    ins_length_dict = {}

    genome_ref = opts.ref_genome_path
    genome_qry = opts.qry_genome_path
    file_path = opts.genome_variation_folder
    alignment_file = align_ref_qry(genome_ref, genome_qry)
    querys, refs = read_genome_align_res(alignment_file)
    for query, ref in zip(querys, refs):
        cal_variation(query, ref, sub_base_dict, sub_length_dict, del_length_dict, ins_base_dict, ins_length_dict)
    write_genome_variation(sub_base_dict, sub_length_dict, del_length_dict, ins_base_dict, ins_length_dict, file_path)

if __name__ == '__main__':
    main()
