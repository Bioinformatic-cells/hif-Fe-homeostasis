#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 14 11:08:20 2019

@author: piazzai
"""
import time

import os  # to manipulate directory

import shutil  # to copy files

import pandas as pd

import numpy as np

import config as config

def log_file(
        run_path,
        script,
        vfile,
        pfile,
        rfile,
        mode,
        nb_pts,
        FORCED_SAMPLING_SCALE,
        RESAMPLING_SIZE,
        skip,
        magnitude,
        increment,
        MAX_ITERATION_SS_IMPROVEMENT,
        SIMU_RTOL,
        SIMU_ATOL,
        EPSILON,
        NB_STEADY_POINTS,
        NUMERICAL_SOLVER,
        Name_Model,
        VAR_NAMES,
        PARAM_NAMES,
        PARAM_ALGEB_DEFS,
        NB_DIM_TOT,
        NB_SAMPLED_DIM,
        DEFAULT_ADD_PARAMS,
        NB_VARS,
        MULTIPLE_INPUTS,
        Fe_GROWTH,
        VARS_BOUNDS,
        PARAMS_CONSTRAINT_BOX,
        Fe_By_Ft_CONSTRAINT,
        MAX_Fe_CONCENTRATION,
        CSTRT_Fetot_HIGH_END,
        CSTRT_Fetot_LOW_END,
        MIN_FREE_IRON,
        NB_RATIOS,
        START,
        INC_m,
        END_SS_m,
        END_MAIN_m,
        constraints_on_Fe_sto,
        constraints_on_TfR1,
        constraints_on_Ft,
        Constraints_on_Ft0,
        Constraints_on_Fetot,
        upperbound_on_C,
        list_sols,
        index_sol,
        list_solutions):


    # IP 13 03 19 begin
    # ================= create files and file info ========================
    shutil.copy(vfile, run_path)
    shutil.copy(pfile, run_path)
    new_info_file = os.path.join(run_path, 'info.txt')
    ifile = open(new_info_file, 'a')
    ifile.write('Execution information of ' +
                time.strftime("%d_%m_%Y-%Hh%Mm%Ss") + '\n')
    ifile.write('Original inputs for variables bounds file' +
                os.path.realpath(vfile) + '\n')
    ifile.write('Original inputs for parameters bounds file' +
                os.path.realpath(pfile) + '\n')
    ifile.write('Script name: ' + script + '\n')
    ifile.write('Input variables file: ' + vfile + '\n')
    ifile.write('Input parameters file: ' + pfile + '\n')
    ifile.write('output results file path: ' + rfile + '\n')
    ifile.write('-----------------------------------------------------' +
                'Settings ' +
                '-----------------------------------------------------\n')
    ifile.write('Mode: ' + mode + '\n')
    ifile.write('Number of iterations: ' + str(nb_pts) + '\n')
    ifile.write('---------------------------------------------' +
                'Steady state information ' +
                '---------------------------------------------\n')
    ifile.write('The steady state was observed during the time window ' +
                'of: ' + str(MAX_ITERATION_SS_IMPROVEMENT) +
                ' iterations \n')
    ifile.write('The relative tolerance of the numerical solver is: ' +
                str(SIMU_RTOL) + '\n')
    ifile.write('The absolute tolerance of the numerical solver is: ' +
                str(SIMU_ATOL) + '\n')
    ifile.write('the maximum authorized variation of a "steady" ' +
                'variable used in check_stability is EPSILON: ' +
                str(EPSILON) + '\n')
    ifile.write('Number of point with quasi-constant value for consider ' +
                'the quasi-steady state is NB_STEADY_POINTS: ' +
                str(NB_STEADY_POINTS) + '\n')
    ifile.write('The numerical solver used is: ' + NUMERICAL_SOLVER +
                '\n')
    ifile.write('The scale sampling is: ' + FORCED_SAMPLING_SCALE + '\n')
    ifile.write('The resampling size is: ' + str(RESAMPLING_SIZE) + '\n')
    ifile.write('The skipped points for Sobol sequence are: ' + str(skip) + '\n')
    ifile.write('In Shell-Sampling mode magnitudine is: ' + str(magnitude) + '\n')
    ifile.write('In Shell-Sampling mode increment is: ' + str(increment) + '\n')
    ifile.write('---------------------------------------------' +
                'Model ' +
                '---------------------------------------------\n')
    ifile.write('The model used is: ' + Name_Model + '\n')
    ifile.write('The variable for that model are ' + str(NB_VARS) + ':\n')

    ifile.write('=====================\n')
    ifile.write(pd.DataFrame(VAR_NAMES, np.arange(1, 8), columns=[
            "Varieble name"]).to_string() + '\n')
    ifile.write('=====================\n\n')

    ifile.write('The parameters for that model are:\n')
    ifile.write('=====================\n')
    ifile.write(pd.DataFrame(PARAM_NAMES, np.arange(1, len(PARAM_NAMES)+1), columns=[
            "Varieble name"]).to_string() + '\n')
    ifile.write('=====================\n\n')

    ifile.write('The parameters algebrically defined are:\n')
    ifile.write('=====================\n')
    ifile.write(pd.DataFrame(PARAM_ALGEB_DEFS, np.arange(1, 4), columns=[
            "Varieble name"]).to_string() + '\n')
    ifile.write('=====================\n\n')
    ifile.write('The total dimension of the parameters space (sampled ' +
                'and not sampled parameters) is: ' + str(NB_DIM_TOT) +
                '\n')
    ifile.write('Dimension of the sampled parameters space is: ' +
                str(NB_SAMPLED_DIM) + '\n')
    if DEFAULT_ADD_PARAMS == [[]]:
        ifile.write('There are no fixed parameters')
    else:
        ifile.write('Default values list of fixed parameters is: ' +
                    str(DEFAULT_ADD_PARAMS) + '\n')
    ifile.write('The concentrations values used for plot the time ' +
                'evolution of our variables are at:\n')
    # MULTIPLE_INPUTS = '{:.2f}'.format(MULTIPLE_INPUTS)
    ifile.write('=================================\n')
    ifile.write(pd.DataFrame(MULTIPLE_INPUTS, 
                             np.arange(1, len(config.MULTIPLE_INPUTS)+1),
                             columns=[
                    "Concentration value (nano Molarity)"]).to_string() +
                '\n')
    ifile.write('=================================\n\n')
    ifile.write('The value of the external iron concentration in the' +
                'Growing case is: ' + str(Fe_GROWTH) + '\n')
    ifile.write('The constraint variales bounds on the initial' +
                'conditions are:\n')

    ifile.write('==========================================\n')
    ifile.write(pd.DataFrame(
            VARS_BOUNDS, VAR_NAMES, columns=[
                    "low bound", "Up bound"]).to_string() + '\n')
    ifile.write('==========================================\n\n')
    ifile.write('The constraint parameters bounds on the initial' +
                'conditions are:\n')
    ifile.write('===========================================\n')
    ifile.write(pd.DataFrame(
            PARAMS_CONSTRAINT_BOX, np.arange(1, len(
                    PARAMS_CONSTRAINT_BOX)+1), columns=[
                    "low bound", "Up bound"]).to_string() + '\n')
    ifile.write('===========================================\n\n')
    ifile.write('The maximum iron ammount load by Ferritin is: ' +
                str(Fe_By_Ft_CONSTRAINT) + '\n')
    ifile.write('The maximum iron concentration during the whole' +
                'cycle is: ' + str(MAX_Fe_CONCENTRATION) + '\n')
    ifile.write('Upper bounds on the [Total Iron-residual Iron]' +
                'at the end: ' + str(CSTRT_Fetot_HIGH_END) + '\n')
    ifile.write('Lower bounds on the [Total Iron-residual Iron]' +
                'at the end: ' + str(CSTRT_Fetot_LOW_END) + '\n')
    ifile.write('Lower bounds on the Iron at the end of iron' +
                'deprivation: ' + str(MIN_FREE_IRON) + '\n')
    ifile.write('Number of ratios for e for each experiment): ' +
               str(NB_RATIOS) + '\n')
    ifile.write('---------------------------------------------' +
                'Setting of time evolution ' +
                '---------------------------------------------\n')
    ifile.write('Initial time for the main simulation ' +
                'at: ' + str(START) + ' hours\n')
    ifile.write('Time step increment of the simulation is: ' +
                str(INC_m) + ' minutes\n')
    ifile.write('Final time for the steady state search is ' +
                'at: ' + str(END_SS_m) + ' hours\n')
    ifile.write('Final time for the main simulation ' +
                'at: ' + str(END_MAIN_m) + ' hours\n')
    ifile.write('---------------------------------------------' +
                'Setting of filter ' +
                '---------------------------------------------\n')
    ifile.write('Fe-sto constrain are: ' '\n')
    ifile.write('==========================================\n')
    ifile.write(pd.DataFrame(
            constraints_on_Fe_sto, [
                    'ratio1', 'ratio2', 'ratio3'], columns=[
                    "low bound", "Up bound"]).to_string() + '\n')
    ifile.write('==========================================\n\n')

    ifile.write('TfR1 constrain are: ' '\n')
    ifile.write('==========================================\n')
    ifile.write(pd.DataFrame(
            constraints_on_TfR1, ['ratio1', 'ratio2', 'ratio3'], columns=[
                    "low bound", "Up bound"]).to_string() + '\n')
    ifile.write('==========================================\n\n')

    ifile.write('Ft constrain are: ' '\n')
    ifile.write('==========================================\n')
    ifile.write(pd.DataFrame(
            constraints_on_Ft, ['ratio1', 'ratio2', 'ratio3'], columns=[
                    "low bound", "Up bound"]).to_string() + '\n')
    ifile.write('==========================================\n\n')
    ifile.write('Ft0 constrains are: ' + str(Constraints_on_Ft0) + '\n')
    ifile.write('Fetot constrains are: ' + str(Constraints_on_Fetot) +
                '\n')
    ifile.write('C upper limit is: ' + str(upperbound_on_C) + '\n')
    ifile.write('==========================================\n\n')
    ifile.write('simulate found: ' + str(len(list_sols)) + '\n')
    i = 0
    if len(list_sols) >= 1:
        list_par = [col[i] for col in list_sols]
        list_parT = list(zip(*list_par))
        name_col = ['set_' + str(jj) for jj in range(1,
                    len(list_sols) + 1)]
        list_parT = pd.DataFrame(
                list_parT, index=config.PARAM_NAMES, columns=name_col)
        list_parT.values.sort()
        tab_par_b_sols = list(zip(
                list_parT['set_1'],
                list_parT['set_' + str(len(list_sols))],
                config.PARAM_NAMES))
        ifile.write(pd.DataFrame(
                tab_par_b_sols, np.arange(1, config.NB_DIM_TOT+1),
                columns=['lower bounds', 'upper bounds',
                         'name']).to_string())
    sol = list(zip(*list_solutions))
    sol = np.array(sol)
    sol_filtered = sol[index_sol]
    if len(sol_filtered) >= 1:
        sol_filteredT = list(zip(*sol_filtered))
        name_col = ['set_' + str(jj) for jj in range(1, len(sol_filtered) + 1)]
        list_par_filter_T = pd.DataFrame(sol_filteredT,
                                         index=config.PARAM_NAMES,
                                         columns=name_col)
        list_par_filter_T.values.sort()
        tab_par_b_solsF = list(zip(
                list_par_filter_T['set_1'],
                list_par_filter_T['set_' + str(len(sol_filtered))],
                config.PARAM_NAMES))
        ifile.write('\n==========================================\n\n')
        ifile.write('Filter found: ' + str(len(sol_filtered)) + '\n')
        ifile.write(pd.DataFrame(tab_par_b_solsF,
                                 np.arange(1, config.NB_DIM_TOT+1),
                                 columns=['lower bounds', 'upper bounds',
                                          'name']).to_string())
    ifile.close()
