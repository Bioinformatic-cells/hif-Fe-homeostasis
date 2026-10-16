#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec  9 16:37:05 2019

@author: piazzai
"""

#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Dec  9 14:56:23 2019

@author: piazzai
"""
# Import the necessary packages and modules
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import config as config
import sys
import inspect
import os  # to manipulate directory
import time

# path information and creation of working path per run
path = os.getcwd()
print('The current working directory is %s' % path)
run_path = path + "/Plot_correlation_params_" + time.strftime("%d_%m_%Y-%Hh%Mm%Ss") + "/"
print(run_path)
# path = "/home/piazzai/Datas/Alavoro_Post-doc1/Simulate/Experiment"

try:
    os.mkdir(run_path)
except OSError:
    print(("Creation of the directory %s failed" % run_path))
else:
    print(("Successfully created the directory %s " % run_path))

path_plot = run_path+'/plot_correlation_params'

try:
    os.mkdir(path_plot)
except OSError:
    print(("Creation of the directory %s failed" % path))
else:
    print(("Successfully created the directory %s " % path))

# _____________________________________________________________________________
filename1 = sys.argv[1]
filename2 = sys.argv[2]

def main():
    plot_correlation_param2(filename1,filename2)
#    print params_listT, '\n'


def plot_correlation_param2(filename1, filename2):
    try:
        data1 = pd.read_csv(filename1,
                           header=0,
                           index_col=0,
                           delim_whitespace=True)
        data2 = pd.read_csv(filename2,
                           header=0,
                           index_col=0,
                           delim_whitespace=True)
        print(data1)
        print(data2)
    except IOError:
        func_name = inspect.currentframe().f_code.co_name
        print('\n ERROR in %s: cannot open file %s\nabort...'\
            % (func_name, filename1), file=sys.stderr)
        sys.exit(1)
    data_paramT1 = data1[0:config.NB_DIM_TOT].values
    data_paramT2 = data2[0:config.NB_DIM_TOT].values
    N_dim = len(data_paramT1)
#    data_paramT = zip(*data_p)
    for ii in range(0,N_dim-1):
        for jj in range(ii+1, N_dim):

            plt.scatter(data_paramT2[ii],data_paramT2[jj], c='r', marker='.', label='Filter_sol_25_02_2020-01h18m28s')
            plt.scatter(data_paramT1[ii],data_paramT1[jj], c='b', marker='.', label='Filter_sol_28_02_2020-04h33m50s')
            plt.legend(loc='upper left')
#            plt.xlim([min(data_paramT[ii])-1*min(data_paramT[ii]),
#                      max(data_paramT[ii])+0.05*max(data_paramT[ii])])
#            plt.ylim([min(data_paramT[jj])-1*min(data_paramT[jj]),
#                       max(data_paramT[jj])+0.05*max(data_paramT[jj])])
            plt.xlim([min(data_paramT2[ii])*0.1, max(data_paramT2[ii])*10])
            plt.ylim([min(data_paramT2[jj])*0.1,max(data_paramT2[jj])*10])
            plt.xlim([min(data_paramT1[ii])*0.1, max(data_paramT1[ii])*10])
            plt.ylim([min(data_paramT1[jj])*0.1,max(data_paramT1[jj])*10])
            plt.xlabel(config.PARAM_NAMES[ii])
            plt.ylabel(config.PARAM_NAMES[jj])
            plt.yscale('log')
            plt.xscale('log')
            plt.gcf()
            plt.savefig(path_plot + '/' + str(ii)+'-'+str(jj)+'Plot_' +
                        config.PARAM_NAMES[ii] + 'VS' +
                        config.PARAM_NAMES[jj] + '_2SP.png', dpi=200)
            plt.close()
    return



if __name__ == '__main__':
    main()

