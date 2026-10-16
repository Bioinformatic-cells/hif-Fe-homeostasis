#! /usr/bin/env python
# -*- coding: utf-8 -*-

"""
This script takes simulation results from search_simplif.py and filters them
according to the constraints defined below: the ratios for the 3 replenishment
experiments; the initial concentrations of Ft and Fe_tot; the upper bound of C.

At the end the min and max of parameters are computed and written to file.

"""

# The computations relative to the ratio are currently valid in the case of
# THREE replenishment experiments  (and TWO constraints for each experiment,
# the one on Fe_sto being presently switched off).

# Note: This can replace minmaxparams.py  (when no filtering is performed)

# The part after "Is this a solution?" could be put in a function

# Usage:
# python filter.py [simu_results_file]
# The argument is the name of file containing the list of simulation results

import sys
import time
import pandas as pd
import inspect
import os
import matplotlib.pyplot as plt
from config import PARAM_NAMES

from config import constraints_on_Fe_sto, constraints_on_TfR1
from config import constraints_on_Ft, Constraints_on_Ft0, Constraints_on_Fetot
from config import upperbound_on_C, NB_DIM_TOT, NB_SAMPLED_DIM

# Name of output files
# IP 14 03 19 I changed the name of the main function in order to be calleb by
# simulate.py si
#     Get name of input file
filename = None
run_path = None


def main():
    global filename, run_path
    filename = sys.argv[1]
    run_path = sys.argv[2]
    run_path, path_plot_correlation_param = directory()
    save_index_sol, list_sols = Filter(filename, run_path)
    plot_correlation_param_filter_here(list_sols, save_index_sol,path_plot_correlation_param)
    print(save_index_sol)

def Filter(filename, run_path):

    print('\n# Simulation results from', filename, '\n')
    # IP name of output files
    param_box_filename = run_path + '/Filter_box_' + time.strftime(
            "%d_%m_%Y-%Hh%Mm%Ss")+'.txt'
    param_box_filename2 = run_path + '/Filter_box_for_input_' + time.strftime(
            "%d_%m_%Y-%Hh%Mm%Ss")+'.txt'
    param_list_filename = run_path + '/Filter_sol_' + time.strftime(
            "%d_%m_%Y-%Hh%Mm%Ss")+'.txt'

    # Open and read the input file
    with open(filename) as f:
        first_line = f.readline()
        if first_line == 'set_1':
            try:
                data = pd.read_csv(filename,
                                   header=0,
                                   index_col=0)

            except IOError:
                func_name = inspect.currentframe().f_code.co_name
                print('\n ERROR in %s: cannot open file %s\nabort...'\
                    % (func_name, filename), file=sys.stderr)
                sys.exit(1)
        else:
            try:
                data = pd.read_csv(filename,
                                   header=0,
                                   index_col=0,
                                   delim_whitespace=True)

        #        print data
            except IOError:
                func_name = inspect.currentframe().f_code.co_name
                print('\n ERROR in %s: cannot open file %s\nabort...'\
                    % (func_name, filename), file=sys.stderr)
                sys.exit(1)
#    with open(filename, 'r') as file:
    Nsols = 0
    newcol = []
    save_index_sol = []
    for i in range(1,len(data.values[0,:])):
        sat = 1
#        print data['set_'+str(i)]
        s = check_ratios(0,
                         data['set_'+str(i)]['ratio_Fe-sto_250'],
                         data['set_'+str(i)]['ratio_TfR1_250'],
                         data['set_'+str(i)]['ratio_Ft_250'],
                         constraints_on_Fe_sto,
                         constraints_on_TfR1,
                         constraints_on_Ft)
#       sat = 0 if one of them is unsat
        if s == 0:
            sat = 0
            continue

        s = check_ratios(1,
                         data['set_'+str(i)]['ratio_Fe-sto_25'],
                         data['set_'+str(i)]['ratio_TfR1_25'],
                         data['set_'+str(i)]['ratio_Ft_25'],
                         constraints_on_Fe_sto,
                         constraints_on_TfR1,
                         constraints_on_Ft)
#       sat = 0 if one of them is unsat
        if s == 0:
            sat = 0
            continue
        s = check_ratios(2,
                         data['set_'+str(i)]['ratio_Fe-sto_10'],
                         data['set_'+str(i)]['ratio_TfR1_10'],
                         data['set_'+str(i)]['ratio_Ft_10'],
                         constraints_on_Fe_sto,
                         constraints_on_TfR1,
                         constraints_on_Ft)
#       sat = 0 if one of them is unsat
        if s == 0:
            sat = 0

      # sat = 0 if one of them is unsat
        if data['set_'+str(i)]['Ft'] < Constraints_on_Ft0[0] or data[
            'set_'+str(i)]['Ft'] > Constraints_on_Ft0[1]:
            sat = 0
        if data['set_'+str(i)]['Fe_Tot_250'] < Constraints_on_Fetot[1] or data[
            'set_'+str(i)]['Fe_Tot_250'] > Constraints_on_Fetot[2]:
            sat = 0
        if data['set_'+str(i)]['Fe_Tot_0'] > Constraints_on_Fetot[0]:
            sat = 0
        if data['set_'+str(i)]['C'] > upperbound_on_C[0]:
            sat = 0
        # ----- Is this a solution?
        if sat == 1:
            Nsols += 1
            save_index_sol.append(i-1)
            newcol.append(data['set_'+str(i)].values.tolist())
#            print newcol, '\n'

    if Nsols >= 1:
        print('Filter returned a number of solution =', Nsols)
        indice = data.index.values.tolist()
        name_col = ['set_' + str(jj) for jj in range(1, Nsols + 1)]
        newcolT = list(zip(*newcol))
        with open(param_list_filename, 'a') as fileR:
            fileR.write(pd.DataFrame(
                newcolT, index=indice, columns=name_col).to_string() + '\n')
        filtered_sols_sort = pd.DataFrame(
            newcolT[0:NB_DIM_TOT], index=indice[0:NB_DIM_TOT],
            columns=name_col)
        filtered_sols_sort.values.sort()
        if Nsols > 1:
            filtered_sols_box = list(zip(
                filtered_sols_sort['set_1'],
                filtered_sols_sort['set_' + str(Nsols)]))
            name_col2 = ['lower_bounds', 'upper_bounds']
            col = pd.Index(name_col2, name='name')
            filtered_sols_box = pd.DataFrame(
                filtered_sols_box, index=indice[0:NB_DIM_TOT], columns=col)
            with open(param_box_filename, 'a') as fileRR:
                fileRR.write(filtered_sols_box.to_string() + '\n')
            filtered_sols_box2 = pd.DataFrame(
                filtered_sols_box, index=indice[0:NB_SAMPLED_DIM], columns=col)
            with open(param_box_filename2, 'a') as fileRR2:
                fileRR2.write(filtered_sols_box2.to_string() + '\n')
            print(filtered_sols_box)
        if Nsols == 1:
            print('Is not possible to define a new filtered box')
    else:
        print('\n =========== Filter NO solutions found ============')
    # Keep only parameter rows for correlation plots (drop ratio/FeTot block).
    list_sols = data.iloc[0:NB_DIM_TOT].values
    # Drop the full results table from this scope before returning.
    del data
    return save_index_sol, list_sols

def check_ratios(
        i_replen, r_Fe_sto, r_TfR, r_Ft_tot, constraints_on_Fe_sto,
        constraints_on_TfR1, constraints_on_Ft):
    sat = 1
    # Constraint on Fe_sto  (from fig 3-B of BBA)
    # r_Fe_sto = Fe_sto_52h / Fe_sto_25h
    bound_inf, bound_sup = constraints_on_Fe_sto[i_replen]
    if r_Fe_sto > bound_sup or r_Fe_sto < bound_inf:
        sat = 0

    # Constraints on TfR1 and Ft (from BBA 2015 paper)
    # r_TfR = TfR1_52h / TfR1_25h
    bound_inf, bound_sup = constraints_on_TfR1[i_replen]
    if r_TfR > bound_sup or r_TfR < bound_inf:
        sat = 0

    # r_Ft_tot = Ft_tot_52h / Ft_tot_25h
    bound_inf, bound_sup = constraints_on_Ft[i_replen]
    if r_Ft_tot > bound_sup or r_Ft_tot < bound_inf:
        sat = 0
    return sat

def plot_correlation_param_filter_here(data, saved_index_sol,path_plot_correlation_param):
    N_dim = len(saved_index_sol)

    data_param_filter = data[:,saved_index_sol]
    for ii in range(0,len(data_param_filter)-1):
        for jj in range(ii+1, len(data_param_filter)):
            plt.scatter(data_param_filter[ii],
                        data_param_filter[jj], marker = '.')
            plt.xlim([min(data_param_filter[ii])-
                      0.15*min(data_param_filter[ii]),
                      max(data_param_filter[ii])+
                      0.05*max(data_param_filter[ii])])
            plt.ylim([min(data_param_filter[jj])-
                      0.15*min(data_param_filter[jj]),
                       max(data_param_filter[jj])+
                       0.05*max(data_param_filter[jj])])
            plt.xlabel(PARAM_NAMES[ii])
            plt.ylabel(PARAM_NAMES[jj])
            plt.yscale('log')
            plt.xscale('log')
            plt.gcf()
            plt.savefig(path_plot_correlation_param +'/' + str(ii)+'-'+str(jj)+
                        'Plot_' +
                        PARAM_NAMES[ii] + 'VS' +
                        PARAM_NAMES[jj] + '.png', dpi=200)
            plt.close()
    return

def directory():
        # path information and creation of working path per run
    path = os.getcwd()
    print('The current working directory is %s' % path)
    run_path = path + "/Run_Filter" + time.strftime("%d_%m_%Y-%Hh%Mm%Ss") + "/"
    print(run_path)
    # path = "/home/piazzai/Datas/Alavoro_Post-doc1/Simulate/Experiment"

    try:
        os.mkdir(run_path)
    except OSError:
        print(("Creation of the directory %s failed" % run_path))
    else:
        print(("Successfully created the directory %s " % run_path))

    path_plot_correlation_param = run_path+'/plot_filtered_correlation_params'

    try:
        os.mkdir(path_plot_correlation_param)
    except OSError:
        print(("Creation of the directory %s failed" % path_plot_correlation_param))
    else:
        print(("Successfully created the directory %s " % path_plot_correlation_param))
    return run_path, path_plot_correlation_param

if __name__ == '__main__':
    main()
