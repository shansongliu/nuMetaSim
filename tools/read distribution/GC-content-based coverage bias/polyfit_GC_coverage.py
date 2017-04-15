"""
Function:     get the ployfit of GC content and coverage bias 
"""

import os
import numpy as np
from optparse import OptionParser

def read_GC_coverage(file_name):
    """
    Args:
        file_name:         str; absolute path of "GC_coverage.txt"                
    """
    GC_list = []
    coverage_list = []
    with open(file_name) as f_r:
        f_r.readline()
        for line in f_r:
            line = line.strip()
            GC_coverage = line.split(':')
            for term in GC_coverage:
                GC, coverage = term.split('-')
                if coverage == '0':
                    GC, coverage = np.int(GC), 0
                else:
                    GC, coverage = np.int(GC), np.float(coverage)
                GC_list.append(GC)
                coverage_list.append(coverage)
    return GC_list, coverage_list

def GC_coverage_normalization(GC_list, coverage_list):
    """
    Args:
        GC_list:         list; a list contains the defined GC range (0-100%)
        coverage_list:    list; a list contains the normalized read count corresonding to its GC content            
    """
    GC_range = 101

    # if coverage < 1.0, then force it to 0
    for i in xrange(len(coverage_list)):
        if coverage_list[i] < 1.0:
            coverage_list[i] = 0

    # normalize
    normalizer = np.zeros((np.array(coverage_list).shape[0] / GC_range))
    for i in xrange(len(normalizer)):
        normalizer[i] = np.sum(coverage_list[i * GC_range:(i + 1) * GC_range])
    geo_count = -1
    coverage_list_normalized = []
    for i in xrange(len(coverage_list)):
        if geo_count == len(normalizer):
            break
        if i % GC_range == 0:
            geo_count += 1
        if normalizer[geo_count] != 0:
            coverage_list_normalized.append(coverage_list[i] / normalizer[geo_count])
        else:
            coverage_list_normalized.append(0.0)

    # find the coverage bigger than 0
    index = [ind for ind in xrange(len(coverage_list_normalized)) if coverage_list_normalized[ind] > 0.0]    # abandon 0 coverage
    GC_list_temp = [GC_list[ind] for ind in index]
    coverage_list_temp = [coverage_list_normalized[ind] for ind in index]

    # outlier deletion
    coverage_mean = np.mean(coverage_list_temp)
    coverage_std = np.std(coverage_list_temp)
    upper_limit = coverage_mean + 2 * coverage_std
    lower_limit = coverage_mean - 2 * coverage_std
    index = [ind for ind in xrange(len(coverage_list_temp)) 
             if coverage_list_temp[ind] < upper_limit and coverage_list_temp[ind] > lower_limit]            # abandon the outlier
    GC_list_processed = [GC_list_temp[ind] for ind in index]
    coverage_list_processed = [coverage_list_temp[ind] for ind in index]

    # sort
    GC_coverage_tuple = zip(GC_list_processed, coverage_list_processed)
    GC_coverage_tuple = sorted(GC_coverage_tuple, key = lambda x: x[0])
    for i in xrange(len(GC_coverage_tuple)):
        GC_list_processed[i] = GC_coverage_tuple[i][0]
        coverage_list_processed[i] = GC_coverage_tuple[i][1]

    return GC_list_processed, coverage_list_processed

def GC_coverage_polyfit(GC_list, coverage_list, polyfit_order = 2):
    """
    Args:
        GC_list:         list; a list contains the defined GC range (0-100%)
        coverage_list:    list; a list contains the normalized read count corresonding to its GC content
        polyfit_order:    int; polynomial fitting order, default is 2        
    """
    GC_list, coverage_list = np.array(GC_list), np.array(coverage_list)
    coef = np.polyfit(GC_list, coverage_list, polyfit_order)
    poly1d_coef = np.poly1d(coef)
    GC_list, coverage_list_fit = list(GC_list), list(poly1d_coef(GC_list))
    return GC_list, coverage_list_fit, coef

def GC_coverage_plot(file_name, GC_list, coverage_list, cover_list_fit):
    """
    Args:
        file_name:        str; absolute path of "GC_coverage.txt"
        GC_list:         list; a list contains the defined GC range (0-100%)
        coverage_list:    list; a list contains the normalized read count corresonding to its GC content
        cover_list_fit:    list; a list contains the fitted number of coverage
    """
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages

    plt.style.use('ggplot')
    fig = plt.figure()
    fig.patch.set_alpha(0.5)
    plt = fig.add_subplot(111)
    plt.patch.set_alpha(0.5)
    plt.plot(GC_list, coverage_list, 'o', color = 'b')
    plt.plot(GC_list, cover_list_fit, color = 'r', lw = 2)
    plt.set_xlabel("GC content %")
    plt.set_ylabel("Normalized read count")
    plt.set_title("GC content - Normalized read count")

    file_path = ""
    if '/' in file_name:
        file_path = '/'.join(file_name.split('/')[:-1]) + '/'
    else:
        file_path = "./"
    pdf_save_path = file_path + "GC_coverage.pdf"
    pdf_save = PdfPages(pdf_save_path)
    pdf_save.savefig()
    pdf_save.close()
    print "GC coverage bias relationship plot saved in %s." % pdf_save_path    

def write_polyfit_coef(file_name, polyfit_coef):
    """
    Args:
        file_name:        str; absolute path of "GC_coverage.txt"
        polyfit_coef:    np.array; polynomial fitting coefficient
    """
    file_path = ""
    if '/' in file_name:
        file_path = '/'.join(file_name.split('/')[:-1]) + '/'
    else:
        file_path = "./"
    file_name = file_path + "polyfit_coefficient.txt"
    with open(file_name, 'w') as f_w:
        f_w.write("#polynomial fitting coefficient")
        f_w.write('\n')
        polyfit_coef_list = list(polyfit_coef)
        polyfit_coef = '\t'.join([str(coef) for coef in polyfit_coef_list])
        f_w.write(polyfit_coef)
        f_w.write('\n')
    print "Fitted poly coefficient file saved in %s." % file_name    

def main():
    usageReminder = ("usage: python polyfit_GC_coverage.py [-f] GC_coverage_bias_file [-p] polyfit_order")
    parser = OptionParser(usage = usageReminder)
    parser.add_option("-f", action = "store", type = "string", dest = "GC_coverage_bias", \
                      help = "The path of the GC coverage bias file generated by \"cal_GC_coverage_relation.py.\"")        
    parser.add_option("-p", action = "store", type = "string", dest = "polyfit_order", \
                      help = "The polyfit order for fitting a curve. Reconmmended and default order is 2.")
    (opts, args) = parser.parse_args()
    if opts.polyfit_order == None:
        opts.polyfit_order = 2

    GC_coverage_path = opts.GC_coverage_bias
    polyfit_order = opts.polyfit_order

    GC_list, coverage_list = read_GC_coverage(GC_coverage_path)
    GC_list_processed, coverage_list_processed = GC_coverage_normalization(GC_list, coverage_list)
    GC_list_processed, coverage_list_fit, coef = GC_coverage_polyfit(GC_list_processed, coverage_list_processed, polyfit_order = 2)
    GC_coverage_plot(GC_coverage_path, GC_list_processed, coverage_list_processed, coverage_list_fit)
    write_polyfit_coef(GC_coverage_path, coef)

if __name__ == '__main__':
    main()
