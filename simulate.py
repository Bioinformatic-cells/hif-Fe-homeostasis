#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Module containing all the functions related to numerical simulations.
"""


# TODO: RE-WORK errors handling in these functions
# sample(model.ADD_SAMPLE_FILE,model.NB_ADD_SAMPLES)
# resampling(par_file_name, nb_pts, model.RESAMPLING_SIZE)
# pca_reduced_sampling()

# TODO: Put err_print = err_print.format(..) call on the same line (sometime it
# is only 'err_print.format(..)' )

# TODO: Generate the file names semi-automatically...(using input files name,
# date and number of points for example)

# TODO: Add infos on the simulation context at the beginning of the output
# file....

# TODO: The modes 'Additional-Sampling' and  'PCA-guided-sampling' have not
# been tested (or used) in this new version.
#       (Was tested and worked in previous version)
#        -----> Some typos may remain: try with few points first.

# function os to manipulate path and folder


import sys
import time
import datetime
import warnings
import numpy as np
import scipy
from scipy import integrate, optimize
from sobol_lib import i4_sobol_generate  # i4_sobol
from pca import PCA
import matplotlib.pyplot as plt
from recordclass import recordclass
import inspect
import itertools
import pandas as pd
import os  # to manipulate directory

# from julia import Julia
# from diffeqpy import de
import model as model
import filter as Filter
import read_file as read_file
import config as config
import log
# Import the module [model_name].py containing
# the model specific informations and functions as 'model'

# IP 15 03 19
# I created conf.py where are allocated constant variables
# _______ constant variable of simulate.py ____________________________________
from config import DEBUG_S as DEBUG
# DEBUG = 1
warnings.filterwarnings("error")

# import dill # if you need to save the current session
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------


"""
- MAX_ITERATION_SS_IMPROVEMENT: Maximum number of interation to search for a
steady state
- DEBUG: Debug mode can be:
                             0 -> No traces except the iteration
                             1 -> traces on error and unsatisfiability of
                             constraints
                             2 -> More informations when going through some
                             functions
                             3 -> step by step with possibility to stop the
                             execution in some function instead of continuing
SIMU_RTOL, and SIMU_ATOL   : Relative and Absolute tolerance of the numerical
                             solver. Small tolerance stands for High precision,
                              and low performaces
EPSILON                    : Define the maximum authorized variation of a
                            'steady' variable (usedin check_stability)
NB_STEADY_POINTS           : Define the number of steady points necessary
to consider the system stable
NUMERICAL_SOLVER           : Numerical solver used to solve ODEs, currently
                             only 'scipy-odeint-default' is used (and handled).
FORCED_SAMPLING_SCALE      : Parameter allowing to forced the scale of the
                             sampling whatever the size of the box 'linear' for
                             linear scale sampling only 'log' for log scale
                             sampling only (beware of negative values) '' for
                             no forcing
"""

# CLI / run-folder globals. Populated by setup_run_environment() when the
# script is executed as main. Left empty when the module is imported for tests.
args_list = []
script = ''
vfile = ''
pfile = ''
rfile = ''
path = os.getcwd()
run_path = ''
path_solutions = ''
path_filter_solutions = ''
plot_solutions = ''
plot_filter_solutions = ''
path_plot_correlation_param1 = ''
path_plot_correlation_param = ''
ofolder = ''
g_time = []


def _mkdir(directory):
    try:
        os.mkdir(directory)
    except OSError:
        print(("Creation of the directory %s failed" % directory))
    else:
        print(("Successfully created the directory %s " % directory))


def setup_run_environment(argv=None):
    """Parse CLI args and create the timestamped Run_* working directories."""
    global args_list, script, vfile, pfile, rfile, path, run_path
    global path_solutions, path_filter_solutions, plot_solutions
    global plot_filter_solutions, path_plot_correlation_param1
    global path_plot_correlation_param, ofolder

    args_list = list(sys.argv if argv is None else argv)
    script = args_list[0]
    vfile = args_list[1]
    pfile = args_list[2]
    rfile = args_list[3]

    path = os.getcwd()
    print('The current working directory is %s' % path)
    run_path = path + "/Run_" + time.strftime("%d_%m_%Y-%Hh%Mm%Ss") + "/"
    print(run_path)
    _mkdir(run_path)

    path_solutions = run_path + '/Solutions_y_t'
    _mkdir(path_solutions)
    path_filter_solutions = run_path + '/Solutions_filtered_y_t'
    _mkdir(path_filter_solutions)
    plot_solutions = run_path + '/Plot_solutions_y_t'
    _mkdir(plot_solutions)
    plot_filter_solutions = run_path + '/Plot_filtered_solutions_y_t'
    _mkdir(plot_filter_solutions)
    path_plot_correlation_param1 = run_path + '/plot_correlation_params'
    _mkdir(path_plot_correlation_param1)
    path_plot_correlation_param = run_path + '/plot_filtered_correlation_params'
    _mkdir(path_plot_correlation_param)

    ofolder = run_path
    rfile = ofolder + rfile
    return run_path
# -----------------------------------------------------------------------------


def apply_method(par_file_name, var_file_name, output_file, nb_pts,
                 MODE='Box-Sampling'):

    """
    Apply the method given in MODE.
    Either:
        _Sample and simulate parameters in a box;
        _Simulate a list of solutions while re-sampling a subset of parameters
        previously fixed;
        _Build boxes around previous solutions, sample them, and simulate;
        _Build a PCA oriented box, sample it and simulate, or sample only on
        the principal components and simulate;

    After simulation, the parameter sets which satisfy a set of contraints are
    saved in a file.

    Inputs: String of the file name containing parameters ::par_file_name,
    -> In mode: 'Box-Sampling' this file should contains a box
                representing the parameter space to sample
    -> In mode: 'Additional-Sampling', 'Re-Sampling',
                'PCA-guided-sampling' this file constainsa list of parameter
                sets (refers to example files for format).

    String of the file name containing variable bounds ::var_file_name
    <--No used anymore

    String of the file name that will contain the results ::output_file

Number of points ::nb_pts
    -> In mode: 'Box-Sampling','PCA-guided-sampling', this is the number of
                points to sample.
    -> In mode: 'Additional-Sampling', there is no use to this input
    -> In mode: 'Re-Sampling', this is the number of points to resample in each
                new box

String defining the method to apply ::MODE, = 'Box-Sampling' by default
    -> 'Box-Sampling': Sample the box given in 'par_file_name' and simulate
                       them using the model equations. Simulation results which
                       satify a set of constraints are  kept in the file
                       'output_file'. Some parameters can be fixed to (or a
                       list of ) given values, and are not  sampled: these
                       values are defined in model.py using DEFAULT_ADD_PARAMS'
                       Some parameters can be defined by the value of other
                       parameters using some algebraic equations.
                       The function add_params() in model.py defines this case.
                       The total number of parameters is given in model.py by
                       NB_DIM_TOT. The number of sampled parameters is
                       NB_SAMPLED_DIM given in model.py, and should be equal to
                       the dimension of the box defined in 'par_file_name'
    -> 'Additional-Sampling': Given a list of previous solutions where some
                              parameters were fixed (and not sampled), it now
                              samples these additional parameters, and
                              re-simulates each previous solution for all of
                              these new samples. The number of additional
                              samples are given in model.py with the constant
                              NB_ADD_SAMPLES. Outside of this mode the possible
                              values of these additional parameters can be
                              fixed using the constant 'DEFAULT_ADD_PARAMS'.
                              These parameters are counted in NB_DIM_TOT but
                              are not count in NB_SAMPLED_DIM as they are not
                              sampled by default. The file containing the box
                              to sample for these parameters is given in
                              model.py in ADD_SAMPLE_FILE.
    -> 'Re-Sampling': It encloses each solution point from 'par_file_name' in a
                      box of size 'RESAMPLING_SIZE' defined in model.py. Then
                      it samples 'nb_pts' in each resulting box, and simulates
                      them. Simulation results which satify a set of
                      constraints are  kept in the file 'output_file'.
    -> 'PCA-guided-sampling': Compute the principal components analysis of the
                              set of solutions given in the file
                              'par_file_name'. If 'PCA_REDUCTION' is set to
                              'full' in model.py it computes an oriented box
                              containing the solutions, samples 'nb_pts'
                              points in it, and simulate them. If
                              'PCA_REDUCTION' is set to 'reduced' in model.py
                              it project the solution on the sub-dimensional
                              linear surface defined by the principal
                              eigenvectors. Then 'nb_pts' points are sampled on
                              this surface, and simulated. Simulation results
                              which satify a set of constraints are  kept in
                              the file 'output_file'.

Ouput: None    <==== Need to change when working on main()
"""
# IP 11/03/19 Example to how print documentation using Sphinx
# _____________________________________________________________________________
    apply_method.__doc__
#    print apply_method.__doc__
# _____________________________________________________________________________
    print('The current MODE is ', MODE)
    print('The number of iterations are', nb_pts)
    # ============================== SAMPLING =================================
#  Read the parameter bounds in the file  (each line -> bounds for 1 parameter)
#  Note: the variables are indexed (in var_bounds and param_bounds) according
#  to their order in the input file.
#  Read the concentration bounds
# Currently not used anymore as bound are hard coded in model module
# read_paramlist, ... , pca_reduced_sampling.
    var_bounds, nb_vars, var_names = \
        read_file.read_file(var_file_name, DEBUG, config.map_VAR_NAMES)
    print(var_bounds, '\n')
    print(var_names, '\n')
    print(config.map_VAR_NAMES.keys(), '\n')
    additional_samples = config.DEFAULT_ADD_PARAMS
    if DEBUG > 2:
        print('\nnb_vars =', nb_vars, '\nvar_bounds =', var_bounds,\
            '\nvar_names  =', var_names)
        print('\n Press Return to continue')
        input()

    if MODE == 'Additional-Sampling':  # TODO change mode name ?
        param_list = read_file.read_params_results(par_file_name)
        # TODO
        additional_samples = sample(
                config.ADD_SAMPLE_FILE, config.NB_ADD_SAMPLES, config.skip)

    elif MODE == 'Re-Sampling':  # resampling of new box TODO
        param_list = resampling(par_file_name, nb_pts, config.RESAMPLING_SIZE)

    elif MODE == 'Shell-Sampling':  # resampling of new box TODO

#        param_list = resampling_shell(par_file_name,
#                                      nb_pts, config.RESAMPLING_SIZE)
        shell_sampling(par_file_name, nb_pts, config.RESAMPLING_SIZE)

    elif MODE == 'PCA-guided-sampling':
        param_list = pca_reduced_sampling(
                par_file_name, nb_pts, dimension=config.PCA_REDUCTION)

    elif MODE == 'Box-Sampling':
        scales = []
        # TODO
        param_list = sample(par_file_name, nb_pts, config.skip,
                            sample_scales=scales)
        # TODO
        param_bounds, nb_params, name = \
            read_file.read_file(par_file_name, DEBUG, config.map_PARAM_NAMES)
    elif MODE == 'read-solutions':
        scales = []
        param_list = read_file.read_params_results(par_file_name)
    else:
        print('\nMode: \'', MODE, '\' is not a valid mode.\nPlease use: \
            \'Box-Sampling\', \'PCA-guided-sampling\', \'Re-Sampling\', \
            \'Additional-Sampling\', or no mode at all.')
        print('\n      >>>> Terminating.')
        return

    if DEBUG > 2:
        print('Number of samples: ', len(param_list))
        print('Sample dimension: ', len(param_list[0]))
        print('\n Press Return to continue')
        input()

# *********** PARAMETERS SIMULATIONS OUT LOOP *********************************
    list_sols = []
    # Loop on the parameter points
    it = 0
    N_sol = 0

    for params in param_list:
        # In case of negative parameters, with PCA resampling
        err, err_print = check_negativity(params)
        if err > 0:
            # it +=1*len(additional_samples)  # counting negative points
            if DEBUG > 0:
                print(err_print)
            if DEBUG > 2:
                assert(0)
            continue

        for sub_sample_param in additional_samples:
            it += 1
            sat = True
            if DEBUG == 0:
                print('Iteration: ', it)
            else:
                print('Iteration: ', it)
                # -------- Recuperating Parameters values ------------

            params, err, err_print = convert_to_recordclass(params)
            if err > 0:
                if DEBUG > 0:
                    print('Errors in parameters conversion to recordclass')
                    print(err_print)
                continue

            if DEBUG > 1:
                # IP 26 03 19 I changed the way to print, more readble
                print('=============================================')
                print('Recuperating parameters')
                tab_params = zip(params, config.PARAM_NAMES)
                print(pd.DataFrame(tab_params,
                                   np.arange(1, len(tab_params)+1),
                                   columns=["parameter", "name"]))
                print('=============================================')
            if DEBUG > 2:
                print('\n Press Return to continue')
                input()

            params, err, err_print = add_sub_sample(params, sub_sample_param,
                                                    copy=False)
# If additional sampled parameters which are not sampled by default
# ex: 'FPN1a' related parameters in iron homeostasis model
            # if DEBUG > 2:
            #     print 'add_sub_sample() output params = ', params
            #     print 'len(params) = ', len(params),'\n'
            if err == True:
                if DEBUG > 0:
                    print('Errors while adding additionnal sampled parameters\
                        on top of previous solutions')
                    print(err_print)
                continue

            params, err, err_print = add_params(params, copy=False)
            # If additional not sampled paramete whose values are defined by
            # the sampled ones ex: 'k_IRP_Ft' in homeostasis model
            # if DEBUG > 2:
            #     print 'add_params output params = ', params
            #     print 'len(params) = ', len(params),'\n'
            if err == True:
                if DEBUG > 0:
                    print('Errors while adding not sampled parameters')
                    print(err_print)
                if DEBUG > 2:                       # Anormal behaviour
                    assert(0)
                continue
            # Constraints on parameters
            params, sat, err, err_print = check_params(params, copy=False)
            if err == True:
                if DEBUG > 0:
                    print(err_print)
                if DEBUG > 2:                       # Anormal behaviour
                    assert(0)
                continue
            # IP 12/10/19 remove the constrains.
            if sat == False:
                if DEBUG > 0:
                    print(err_print)
                    print('Parameters do not respect the constraints')
                continue

            if DEBUG > 2:
                print(' check_params output params = ', params)
                print('len(params) = ', len(params), '\n')
                print('\n Press Return to continue')
                input()
            #err,err_print = update_params(params) #TODO keep ?
            init_vars, sat, err, err_print = generate_initial_state(params)
            # Compute steady-state initial conditions using algebraic equation
            # on parameters Also verify and try to improve steadiness

            if err == True:
                if DEBUG > 0:
                    print('Errors in steady state simulation loops')
                    print(err_print)
                if DEBUG > 2:
                    assert(0)
                continue

            if sat == False:
                if DEBUG > 0:
                    print(err_print)
                    print('No steady state found')
                    print('_______________________________________________' + \
                    '__________________________________')
                continue

            sol, sat, err, err_print = simulate_and_check(
                    params, init_vars, config.MULTIPLE_INPUTS)
            # TODO Changed sol Tuple to array to authorize modifs
            # Simulates the model amd check absolute constraints
            # If multiple experiments, or possible input values
            # it performs mulitple simulations and check results
            # ex: Iron re-apport in Iron Homeostasis model
            # No relative constraints check between simulations
            # Output solution in a model dependant format

            # IP 10/10/19 I count the solution e I print them in a file
            if sat == True:
                N_sol = N_sol+1

                if DEBUG > 1:
                    # Optional per-solution dump (debug only); always close file.
                    results_file = os.path.join(
                            run_path, 'results_y'+str(N_sol)+'.txt')
                    name_col = config.VAR_NAMES
                    col = pd.Index(name_col, name='time')
                    with open(results_file, 'a') as fileR:
                        fileR.write(pd.DataFrame(
                            sol[2][0][0],
                            index=sol[2][0][1],
                            columns=col).to_string() + '\n')
                if N_sol == 1:
                    global g_time
                    g_time = np.asarray(sol[2][0][1]).copy()



            if err > 0:
                if DEBUG > 0:
                    print('Errors in main simulation loops')
                    print(err_print)
                if DEBUG > 2:
                    print('CLOSING SIMULATION')
                    assert(0)
                continue

            if sat == False:
                if DEBUG > 0:
                    print('Parameters set is not retained as a potential ' + \
                          'solution')
                    print(err_print)
                continue

            sol[2] = compute_Fe_tot(sol[2])
            # wipe the simulations to free memory when not in 'Plot' mode
            # IP 25/11/2019 and replacing the simulation with the values of
            # Fe_tot

            if err == True:
                if DEBUG > 0:
                    print('Errors saving current solution')
                    print(err_print)
                if DEBUG > 2:
                    assert(0)
                continue
            list_sols.append(sol)

    err, err_print = print_end_infos(list_sols, par_file_name)
    # Print results in file on console
    # Default Output: NB sols and enclosing box
    # If specified in model print additional infos
    # Error handling
    if err == True:
        print('Errors while writing complete simulation sumary ')  # TODO
        if DEBUG > 0:
            print(err_print)
        if DEBUG > 2:
            assert(0)



    return len(list_sols), list_sols
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# TODO error code


def sample(file_name, nb_points, skip, sample_scales=[]):

    """
    Read a file defining the bounding box of a parameter space
    and output a list of sampled points using sobol method

    Inputs: String of the file name::file_name,
            Number of points to be sampled::Npoints

    Ouput: Liste of sampled parameters::param_list

"""
# IP 11/03/19 how print documentation using Sphinx
# _____________________________________________________________________________
    sample.__doc__
#    print sample.__doc__
# _____________________________________________________________________________

    param_bounds, N_par, name = read_file.read_file(file_name,
                                                    DEBUG,
                                                    config.map_PARAM_NAMES)
    params = [0]*N_par
    init = [0]*N_par
    param_list = [init]*nb_points

    contains_neg_bounds, idx = check_negative_bounds(param_bounds)

    # TODO
    # Negative bounds are not an error (exemple sampling on centrered space in
    # pca method, or electric related parameters)
    # However we must have (bu/bl < 100) to avoid the log scale method or
    # update the code for it
    if contains_neg_bounds == False:
        param_bounds_log = [(np.log10(item[0]), np.log10(item[1]), item[2])
                            for item in param_bounds]
    else:
        print('Bound on Parameter ', name[idx[0]], 'is negative: ',\
            param_bounds[idx[0]][idx[1]])
        print('Forcing sampling to be linear')
        # global FORCED_SAMPLING_SCALE
        config.FORCED_SAMPLING_SCALE = 'linear'

        # assert(contains_neg_bounds == False)

    # Draw a list of random points in the unit hyper-cube [0,1]^nb_params
    par = i4_sobol_generate(N_par, nb_points, config.skip)
    # Loop on the parameter points
    for ipt in range(0, nb_points):       # IP 22 03 19
        curr_par = par[0:N_par, ipt]  # current parameters set
        for ipar in range(0, N_par):  # read the initial box file and
            pmin = param_bounds[ipar][0]  # save in pmin and pmax low and up
            pmax = param_bounds[ipar][1]  # bounds

            if (config.FORCED_SAMPLING_SCALE == 'log') or (
                    config.FORCED_SAMPLING_SCALE == '' and
                    pmin > 0 and pmax/pmin > 100.):
                # More than 2 orders of magnitude
                # From [0,1] to log scale (sampling in log scale)
                # print name[ipar],': log scale'  --> #TODO store it
                # sample_scales.append('log')
                # -> use log scale
                interv_i = (param_bounds_log[ipar][1] -
                            param_bounds_log[ipar][0]) * curr_par[ipar] \
                            + param_bounds_log[ipar][0]
                try:
                    params[ipar] = 10**(interv_i)
                    if DEBUG > 2:
                        print("DBG: interv_i =", interv_i, '    \
                        param value =', params[ipar])

                except RuntimeWarning:
                    print("WARNING: overflow encountered in power \
                    10**(interv_i) with interv_i =", interv_i)
                    print("    ipar =", ipar, "  bounds =", param_bounds[ipar])
                    continue  # Skip this parameter point

            else:
                #  From [0,1] to normal scale (sampling in normal scale)
                # print name[ipar],': linear scale'
                # sample_scales.append('linear')
                params[ipar] = (pmax - pmin) * curr_par[ipar] + pmin
                if DEBUG > 2:
                    print("DBG: interv_i =", params[ipar])

            param_list[ipt] = params[:]

    return param_list
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------


def check_negative_bounds(param_bounds):
    """
    Check if there is negative values for the bounds on parameters.
    Used in sample(), to ensure the log can be applied.


    Inputs: Array of pair of bounds::param_bounds

    Ouput: Boolean answering the negativity existence test:do_contain,
           Two dimensional index of the first negative bound::index
"""
# IP 11/03/19 how print documentation using Sphinx
# _____________________________________________________________________________
    check_negative_bounds.__doc__
#    print check_negative_bounds.__doc__
# _____________________________________________________________________________
    i = 0
    index = [-1, -1]
    do_contain = False
    for bounds in param_bounds:
        for j in range(len(bounds)):
            if bounds[j] <= 0:
                do_contain = True
                index = [i, j]
                break
        if do_contain == True:
            break

        i += 1

    return do_contain, index
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------


def resampling(list_points_file, resampling_nb_pts, box_size):
    """
    Resample points around a given list of previous solutions.

    INPUTS:
            File containing the previous solutions         :: list_points_file
            Number of points to resample FOR EACH solution :: resampling_nb_pts
            Size of the box built around a solution        :: box_size
            -> given as a % increase/decrease around the solution

    OUTPUT:
            List of resampled points :: resampled_list
"""
# IP 11/03/19 how print documentation using Sphinx
# _____________________________________________________________________________
    resampling.__doc__
#    print resampling.__doc__
# _____________________________________________________________________________

    list_points = read_file.read_params_results(list_points_file)
    list_points = np.array(list_points)
    if DEBUG > 1:
        print('list_points = \n', list_points)
    list_points = list_points[:, :config.NB_SAMPLED_DIM].tolist()
    sizes = [box_size]*config.NB_SAMPLED_DIM
    resampled_list = []
    temporary_file_name = 'temp-resampling-box.txt'
    for par in list_points:
        box = build_box(par,sizes, temporary_file_name)  #TODO output bl,bu ?
        resampled_list += sample(temporary_file_name,
                                 resampling_nb_pts,
                                 config.skip)

    return resampled_list
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------


def build_box(point, sizes, output_file):
    """
    Give the box over-approximatin a point.
    The surrounding box size is given a percentage x%
    In each dimension pi of p the box is:
                if x < 0.5: [pi-x*pi , pi+x*pi]
                if x > 0.5: 10^[log(pi)-x*log(pi) , log(pi)+x*log(pi)]
                (or opposite if log(pi)<0)
    TODO: Add an order of magnitude size size to define the box

    Inputs: Array defining a point::point,
            Array of the relative sizes (+/- x%) of the box in each
            dimension::sizes,
            file name of the results::output_file


    Output: lower and upper bound of the box: box
# in the considered 'scaling' (normal or log)  \
"""
# IP 11/03/19 how print documentation using Sphinx
# _____________________________________________________________________________
    build_box.__doc__
#    print build_box.__doc__
# _____________________________________________________________________________

    box = []
    for i in range(len(point)):
        v = point[i]

        if sizes[i] < 0.0:
            print('Negative sizes are not allowed', sizes[i])
            assert(sizes[i] >= 0.0)
        elif sizes[i] < 0.5:  # normal scale
            box.append([config.PARAM_NAMES[i], v - sizes[i]*v, v + sizes[i]*v])
        elif sizes[i] >= 0.5:  # log scale
            logv = np.log10(v)
            diffv = sizes[i]*np.log10(v)
            if DEBUG > 0:
                print('point[', i, '] = ', v)
                print('log(point[i]) = ', logv)
                print('sizes[i] = ', sizes[i])
                print('sizes[i]*log(point[i]) = ', diffv)
            if diffv >= 0:
                box.append([
                        config.PARAM_NAMES[i],
                        10**(logv - diffv),
                        10**(logv + diffv)])
            else:
                box.append([
                        config.PARAM_NAMES[i],
                        10**(logv + diffv),
                        10**(logv - diffv)])
        else:
            print(sizes[i])
            assert(0)
    name_col = ['name', 'lower_bounds', 'upper_bounds']
    box_file = 'temp-resampling-box.txt'
    try:
        f_box = open(box_file, 'w')
    except IOError:
        print('\nError in build_box(): cannot open file %s\nabort...' \
                                            % box_file)
        sys.exit(1)
    f_box.write(pd.DataFrame(
            box, columns=name_col).to_string(index=False) + '\n')
    f_box.close()
    tab_box = pd.DataFrame(box, columns=name_col)
    boxV = zip(tab_box['lower_bounds'].values, tab_box['upper_bounds'].values)
    return boxV

# -----------------------------------------------------------------------------


def resampling_shell(list_points_file, resampling_nb_pts, box_size):

    data = read_file.read_params_results(list_points_file)
    list_pointsT = zip(*data)
    list_points_pd = pd.DataFrame(list_pointsT)
    list_points_pd.values.sort()
    box_sort = list_points_pd[0:config.NB_SAMPLED_DIM]
    box_extended = [(box_sort[0].values
                     - config.RESAMPLING_SIZE*box_sort[0])/config.magnitude,
                    (box_sort[len(data)-1].values
                        + config.RESAMPLING_SIZE*box_sort[
                                len(data)-1].values)*config.magnitude]

    box_extendedT = zip(*box_extended)
    box_file = 'box_extended.txt'
    name_col = ['lower_bounds', 'upper_bounds']
    col = pd.Index(name_col, name='name')
    try:
        f_box = open(box_file, 'w')
    except IOError:
        print('\nError in build_box(): cannot open file %s\nabort...' \
                                            % box_file)
        sys.exit(1)
    f_box.write(pd.DataFrame(
            box_extendedT, index=config.PARAM_NAMES[0:config.NB_SAMPLED_DIM],
            columns=col).to_string() + '\n')
    f_box.close()
    if DEBUG > 1:
        print('the extended box is:\n')
        print(box_extendedT)

    box_pd = pd.DataFrame(box_extendedT,
                          index=config.PARAM_NAMES[0:config.NB_SAMPLED_DIM],
                          columns=col)
    print(box_pd)
    par_in_external_box = []
    s=0
    while True:
        config.skip = config.skip + config.increment
        param_list = sample(box_file, config.increment, config.skip)
#       for par in param_list:
        for ii in range(0, len(param_list)):
            sat = False
            par = param_list[ii]
            for jj in range(0, len(par)):
                if par[jj] < list_points_pd[0][jj] or par[
                        jj] > list_points_pd[len(data)-1][jj]:
                    sat = True
                    break
            if sat == True:
                s=s+1
                par_in_external_box.append(par)
                print(str(sys.getsizeof(par_in_external_box)/1024/1024))
                print(s)
                if len(par_in_external_box) == resampling_nb_pts:
                    print(len(par_in_external_box))
                    sys.exit
                    return par_in_external_box
# -----------------------------------------------------------------------------

def shell_sampling(list_points_file, resampling_nb_pts, box_size):

    data = read_file.read_params_results(list_points_file)
    list_pointsT = zip(*data)
    list_points_pd = pd.DataFrame(list_pointsT)
    list_points_pd.values.sort()
    box_sort = list_points_pd[0:config.NB_SAMPLED_DIM]
    b_inf = (box_sort[0].values
                     - config.RESAMPLING_SIZE*box_sort[0])/config.magnitude
    b_sup = (box_sort[len(data)-1].values
                        + config.RESAMPLING_SIZE*box_sort[
                                len(data)-1])*config.magnitude
    box_extended = [b_inf,b_sup]
    D_box = box_sort[len(data)-1].values - box_sort[0].values
    thikness_inf = box_sort[0].values - b_inf.values
    thikness_sup = b_sup.values - box_sort[len(data)-1].values
    PP = box_sort[9].values
    variation_inf = - (PP-box_sort[0].values)- thikness_inf*(PP-box_sort[0].values)/D_box
    variation_sup = (box_sort[len(data)-1].values-PP) + thikness_sup*(PP-box_sort[0].values)/D_box
    all_inf = add_all_comb(PP, variation_inf)
    #print a_all_inf
    all_sup = add_all_comb(PP, variation_sup)
    list_param = all_inf + all_sup
    plot_correlation_param(list_param)
#    list_paramT = zip(*list_param)
#    name_col =['set_' + str(jj) for jj in range(1, len(list_param) + 1)]
#    col = pd.Index(name_col, name='name')
#    file_shell = 'file_shell_sampling.txt'
#    try:
#        file_shell_points = open(file_shell, 'w')
#    except IOError:
#        print '\nError in build_box(): cannot open file %s\nabort...' \
#                                            % file_shell
#        sys.exit(1)
#    file_shell_points.write(pd.DataFrame(
#        list_paramT, index=config.PARAM_NAMES[0:config.NB_SAMPLED_DIM],
#        columns=col).to_string() + '\n')
#    file_shell_points.close()
    box_extendedT = zip(*box_extended)
    box_file = 'box_extended.txt'
    name_col = ['lower_bounds', 'upper_bounds']
    col = pd.Index(name_col, name='name')
    try:
        f_box = open(box_file, 'w')
    except IOError:
        print('\nError in build_box(): cannot open file %s\nabort...' \
                                            % box_file)
        sys.exit(1)
    f_box.write(pd.DataFrame(
            box_extendedT, index=config.PARAM_NAMES[0:config.NB_SAMPLED_DIM],
            columns=col).to_string() + '\n')
    f_box.close()
    if DEBUG > 1:
        print('the extended box is:\n')
        print(box_extendedT)

    box_pd = pd.DataFrame(box_extendedT,
                          index=config.PARAM_NAMES[0:config.NB_SAMPLED_DIM],
                          columns=col)
    print(box_pd)

    return
# -----------------------------------------------------------------------------
def increment_ind(ind, M):       #M is the number of the a elements
    if (ind[len(ind)-1] < M-1):
        ind[len(ind)-1] += 1
        return ind
    else:
        t = increment_ind(ind[0:len(ind)-1],M-1)
        t.append(t[len(t)-1]+1)
        return t

def add_1(a,N, variation):                 #N is the index number
    index = []                             # the number of gruop for make the
    for ii in range(0,N):                  # the combination
        index.append(ii)
    a_all = []
    while (True):
        a_new = a[:]
        for jj in range(0,len(index)):
            a_new[index[jj]] += variation[jj]
        a_all.append(a_new)
        if (index[0] == len(a)-N):
            return a_all
        index = increment_ind(index, len(a))

def add_all_comb(a, variation):
    ret = []
    for ii in range(1, len(a)+1):
        ret= ret + add_1(a,ii,variation)
    return ret

#def resampling_shell3(list_points_file, resampling_nb_pts, box_size):
#
#    data = read_file.read_params_results(list_points_file)
#    list_pointsT = zip(*data)
#    list_points_pd = pd.DataFrame(list_pointsT)
#    list_points_pd.values.sort()
#    box_sort = list_points_pd[0:config.NB_SAMPLED_DIM]
#    box_extended = [(box_sort[0].values
#                    - config.RESAMPLING_SIZE*box_sort[0])/config.magnitude,
#                    (box_sort[len(data) - 1].values
#                        + config.RESAMPLING_SIZE*box_sort[
#                                len(data) - 1].values)*config.magnitude]
#
#    box_extendedT = zip(*box_extended)
#    box_file = 'box_extended.txt'
#    name_col = ['lower_bounds', 'upper_bounds']
#    col = pd.Index(name_col, name='name')
#    try:
#        f_box = open(box_file, 'w')
#    except IOError:
#        print '\nError in build_box(): cannot open file %s\nabort...' \
#                                            % box_file
#        sys.exit(1)
#    f_box.write(pd.DataFrame(
#            box_extendedT, index=config.PARAM_NAMES[0:config.NB_SAMPLED_DIM],
#            columns=col).to_string() + '\n')
#    f_box.close()
#    if DEBUG > 1:
#        print 'the extended box is:\n'
#        print box_extendedT
#
#    box_pd = pd.DataFrame(box_extendedT,
#                          index=config.PARAM_NAMES[0:config.NB_SAMPLED_DIM],
#                          columns=col)
#    print box_pd
#
#    param_list = sample(box_file, resampling_nb_pts, config.skip)
#    # for par in param_list:
#    par_in_external_box = []
#    for ii in range(0, len(param_list)):
#        sat = True
#        par = param_list[ii]
#        for jj in range(0, len(par)):
#            if par[jj] >= list_points_pd[0][jj] and par[jj] <= list_points_pd[
#                    len(data)-1][jj]:
#                sat = False
#                break
#        if sat == True:
#            par_in_external_box.append(par)
#            if len(par_in_external_box) == resampling_nb_pts:
#                print len(par_in_external_box)
#                return par_in_external_box
## -----------------------------------------------------------------------------


def pca_reduced_sampling(list_points_file, resampling_nb_pts,
                         dimension='reduced'):
    """
    Build a PCA oriented box containing previous solutions and sample it, or
    sample only on some principal components. PCA code is in pca.py.
    TODO: put it there ?
    Note: In current version the the sampling Box is increased by a percentage
    given in this function (incr_box_size).

    INPUTS:
            File containing the previous solutions          :: list_points_file
            Number of points to resample                    :: resample_nb_pts
            Method to use, either 'full' dimension or
            'reduced' dimension (default)                   :: dimension

    OUTPUTS:
            list of sampled parameters (in dimension 'config.NB_SAMPLED_DIM'
            even if reduced) :: points_list
"""
# IP 11/03/19 how print documentation using Sphinx
# _____________________________________________________________________________
    pca_reduced_sampling.__doc__
#    print pca_reduced_sampling.__doc__
# _____________________________________________________________________________

    # Given a large set of solutions as input
    # Read inital list of points from file
    initial_list = read_file.read_params_results(list_points_file)
    initial_list = np.array(initial_list)
    # 20% sampling box size increase
    incr_box_size = 0.2

    if DEBUG > 1:
        print('\n\n ##### PCA ORIENTED BOX CONSTRUCTION ###### \n')
        print('NB_SAMPLED_DIM = ', config.NB_SAMPLED_DIM)
        print('initial list = ', initial_list)

    # Interested only in NB_SAMPLED_DIM first parameters
    initial_list = initial_list[:, :config.NB_SAMPLED_DIM]

    # PCA of initial list of points
    # scaled_param are points after variable change
    # variable change (canonic base -> eigenvectors base)
    scaled_points, eigvals, eigvects = PCA(initial_list)

    if dimension == 'reduced':
        # Compute the number of considered principal components
        # Method: Proportion of variance ~ 0.9
        # Dimensionality (number of considered principal components)
        nb_pc = 0
        # Current sum of eigenvalues
        cseig = 0
        # Total sum of eigenvalues
        seig = np.sum(eigvals)
        for e in eigvals:
            # Stop when we accumulated 90% of the variance
            if cseig/seig > 0.9:
                break
            cseig += e
            nb_pc += 1
        if DEBUG > 0:
            print('Proportion of Variance = ', cseig/seig)

    elif dimension == 'full':
        if DEBUG > 0:
            # Still compute Pov for debug infos
            # Dimensionality (number of considered principal components)
            # Current sum of eigenvalues
            nb_pc = 0
            cseig = 0
            # Total sum of eigenvalues
            seig = np.sum(eigvals)
            for e in eigvals:
                # Stop when we accumulated 90% of the variance
                if cseig/seig > 0.9:
                    break
                cseig += e
                nb_pc += 1
            print('Proportion of Variance = ', cseig/seig)
            print('Number of principal components = ', nb_pc)
            print('System forced on all dims (', config.NB_SAMPLED_DIM, ')')
        nb_pc = config.NB_SAMPLED_DIM

    else:
        print(dimension, ' is an incorrect argument,\
        choose either \'full\' or \'reduced\'')
        assert(0) #TODO error code

    if DEBUG > 1:
        print('nb_pc = ', nb_pc)

    # Building Orthogonal projection matrix
    # on the nb_pc first directions
    P = np.identity(len(eigvals))
    for i in range(len(eigvals)):
        if i >= nb_pc:
            P[i, i] = 0.0

    # Projections of the points
    proj_points = np.dot(scaled_points, P)

    temporary_file_name = 'pca_sample.box'
    # Ecriture de la boite (non-zeros) dans un fichier temp
    list2box(proj_points[:, :nb_pc], temporary_file_name, incr_box_size)
    # resampling sous dimension
    trans_points_list = sample(temporary_file_name, resampling_nb_pts,
                               config.skip)
    trans_points_list = np.array(trans_points_list)

    # rajout zeros pour arriver a la dimenson probleme
    temp = (trans_points_list.T).tolist()
    nb_zero_par = config.NB_SAMPLED_DIM - nb_pc
    addentum = [[0.0]*resampling_nb_pts] * nb_zero_par
    # Add at the end (eigvects have been ordered)
    temp = temp + addentum
    trans_points_list = (np.array(temp)).T
    ###

    # Inverse of the variable change
    # (eigenvectors base -> canonic base)
    points_list = np.dot(eigvects, trans_points_list.T).T
    points_list *= np.std(initial_list, 0)  # denormalization
    points_list += np.mean(initial_list, 0)  # decentrage
    points_list = points_list.tolist()

    return points_list
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------


def list2box(liste, output_file, size):
    """
   Write in a file the box enclosing a list of paramter sets.
   This box can be bloated using the argument 'size'.

   Inputs: list of parameter sets               ::liste,
           file_name                            ::output_file
           box size increase (percentage)       :: size
   Output: tuple (box_lower_bounds,box_upper_bounds):: (bl,bu)
"""
# IP 11/03/19 how print documentation using Sphinx
# _____________________________________________________________________________
    list2box.__doc__
#    print list2box.__doc__
# _____________________________________________________________________________

    dim = len(liste[0])
    if DEBUG > 0:
        print('DIMENSION = ', dim)

    npl = np.array(liste)
    bl = [np.amin(npl[:, i]) for i in range(dim)]
    bl = [bl_i-size*bl_i for bl_i in bl]
    bu = [np.amax(npl[:, i]) for i in range(dim)]
    bu = [bu_i-size*bu_i for bu_i in bu]

    try:
        f_box = open(output_file, 'w')
    except IOError:
        print('\nError in list2box(): cannot open file %s\nabort...' \
                                            % output_file)
        sys.exit(1)

    longest_string = 0
    for item in config.PARAM_NAMES:
        if len(item) > longest_string:
            longest_string = len(item)

    for i in range(len(bl)):
        if DEBUG > 0:
            print('{vname: <{width}} :  {min}    {max}'\
                .format(vname=config.PARAM_NAMES[i],
                        width=longest_string,
                        min=bl[i],
                        max=bu[i]))
        f_box.write('{vname: <{width}} :  {min}    {max}\n'
                    .format(vname=config.PARAM_NAMES[i],
                            width=longest_string,
                            min=bl[i],
                            max=bu[i]))

    f_box.close()
    return (bl, bu)
# -----------------------------------------------------------------------------

# ######### TODO  ######## ABOVE FUNCTIONS ARE NOT UP TO STANDARD ERROR OUTPUT

# -----------------------------------------------------------------------------
# TODO return 'sat' in addition to 'err'


def check_negativity(p):
    """
    Check if a point p is negative

    INPUT:
            the point :: p

    OUTPUT
            Boolean checking the negativity     :: err
            Some informations in string format  :: err_print
"""
# IP 11/03/19 how print documentation using Sphinx
# _____________________________________________________________________________
    check_negativity.__doc__
#    print check_negativity.__doc__
# _____________________________________________________________________________

    # notice the use of 'err' instead of 'sat'
    err = 0
    err_print = 'No negative value in current point'

    idx = -1
    pi = -1.0
    if DEBUG > 1:
        # IP 26 03 19 I changed the way to print, more readble
        tab_p = zip(p, config.PARAM_NAMES)
        print('____________________________________________')
        print('current point in negativity check:')
        print(pd.DataFrame(tab_p,
                           np.arange(1, len(tab_p)+1),
                           columns=["parameter", "name"]))
        print('____________________________________________')
    for i in range(len(p)):
        if p[i] <= 0:
            err = 1
            idx = i
            pi = p[i]
            break
    if err == 1:
        err_print = 'Point with negative value in {}th direction: p[{}] = {}'
        err_print = err_print.format(idx, idx, pi)

    return err, err_print
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# TODO Handling possible errors ?


def convert_to_recordclass(params):
    # IP 10/12/19 I removed input. I pass the multiple inputs as args of
    # diff_eqs(), simulate(), integrate_scipy_default() and
    # scipy.integrate.odeint ()
    err = 0
    err_print = 'No error in convert_to_recordclass'
    names = ''
    for name in config.PARAM_NAMES:
        names += name + ' '
    named_params = recordclass('p', names)
    params = [params[i] if i < config.NB_SAMPLED_DIM
              else 0.0 for i in range(config.NB_DIM_TOT)]
    params = named_params(*params)

    return params, err, err_print

# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------


def add_sub_sample(params, added_params, copy=True):
    """
    Define the value of the not sampled parameters with a value fixed in
    'DEFAULT_ADD_PARAMS'. These parameters are written after the
    'config.NB_SAMPLED_DIM' sampled parameters.

    INPUTS:
            The parameter set before adding the additional fixed
            paramters                                           :: params
            The value of the added parameters                   :: added_params
            Boolean for in-place or copy                        :: copy

    OUTPUT:
            The completed parameter set     :: params
            The boolean checking the error  :: err
            The associated print            :: err_print
"""
# IP 11/03/19 how print documentation using Sphinx
# _____________________________________________________________________________
    add_sub_sample.__doc__
#    print add_sub_sample.__doc__
# _____________________________________________________________________________

    if copy == True:
        params = params.copy()

    err_print = 'No error'
    err = 0

    size_added_params = len(added_params)
    nb_total_sampled = config.NB_SAMPLED_DIM + size_added_params

    if len(params) == config.NB_SAMPLED_DIM:     # TODO remove if stable
        err = 1
        err_print = 'The number of parameters should be '
        + str(config.NB_DIM_TOT) + ' not ' + str(len(params))
        # params = params + added_params
        # Params is the incomplete list of the default sampled parameters
        # Rq: '+' here is list concatenation
        if DEBUG > 2:  # Anormal behaviour
            print(err_print)
            assert(0)
        return params, err, err_print
    # Should be the normal number of parameteres at the point of the method
    # (21/09/2018)
    elif len(params) == config.NB_DIM_TOT:
        if config.NB_DIM_TOT < nb_total_sampled:
            err = 1
            err_print = 'Total number of parameters \
                in model too low.\n NB_DIM_TOT =  {} < {} = {}+{}'
            err_print = err_print.format(config.NB_DIM_TOT,
                                         nb_total_sampled,
                                         config.NB_SAMPLED_DIM,
                                         size_added_params)
            if DEBUG > 2:  # Anormal behaviour
                print(err_print)
                assert(0)
        elif size_added_params == 0:
            # Nothing to inject (DEFAULT_ADD_PARAMS = [[]])
            pass
        else:
            try:
                # Params is already the full list of parameters; update values
                for offset, value in enumerate(added_params):
                    params[config.NB_SAMPLED_DIM + offset] = value
            except Exception:
                err = 1
                err_print = 'Error in array assignement'
                if DEBUG > 2:  # Anormal behaviour
                    print(err_print)
                    assert(0)
    else:
        err = 1
        err_print = 'This number of input parameters \
            is anormal.\n >>> Use {} or {}'
        err_print = err_print.format(config.NB_SAMPLED_DIM, config.NB_DIM_TOT)
        if DEBUG > 2:  # Anormal behaviour
            print(err_print)
            assert(0)

    return params, err, err_print
# -----------------------------------------------------------------------------

#  -----------------------------------------------------------------------------


def add_params(params, copy=True):
    """
    Define the value of the parameters determined by algebraic equations on
    the other parameters.
    Actual definition is given in the model.add_params() function

    INPUTS:
            Incomplete parameter set :: params
            In-place or copy mode    :: copy

    OUTPUT:
            Completed parameter set  :: params
            Error checking boolean   :: err
            Associated print         :: err_print
"""
# IP 11/03/19 how print documentation using Sphinx
# _____________________________________________________________________________
    add_params.__doc__
#    print add_params.__doc__
# _____________________________________________________________________________

    params, err, err_print = model.add_params(params, copy)
    return params, err, err_print
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------


def check_params(params, copy=True):
    """
    Check some constraints on the parameters. Can modify the parameters if
    needed.
    Actual model specific code in  model.check_static_constraints_params()

    INPUTS:
            Parameter set to be checked     :: params
            In-place or copy mode           :: copy

    OUTPUT
            Parameter set if modified               :: params
            Boolean check if constraints passed     :: sat
            Boolean check if errors                 :: err
            Print associated to errors and/or unsat :: err_print
"""
# IP 11/03/19 how print documentation using Sphinx
# _____________________________________________________________________________
    check_params.__doc__
#    print check_params.__doc__
# _____________________________________________________________________________

    if copy == True:
        params = params.copy()
    params, sat, err, err_print = \
        model.check_static_constraints_params(params, copy)
    return params, sat, err, err_print
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------


def update_params(params):  # TODO keep ? (Not used currently)
    """
    Prototype of function to pass the parameters as global variable in the
    model module
"""
# IP 11/03/19 how print documentation using Sphinx
# _____________________________________________________________________________
    update_params.__doc__
#    print update_params.__doc__
# _____________________________________________________________________________

    err, err_print = model.set_params(params)
    return err, err_print
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------


def generate_initial_state(params):
    """
    Given a parameter set, it generates the associated initial conditions.
    Then it check if the initial conditions are a steady state, and try to
    simulate and reach one if not.
    Then static constraints are verified on the obtained steady state.

    INPUT:
            Parameter set   :: params

    OUTPUTS:
            Initial conditions                                  :: init_vars
            Boolean checking validity of the initial conditions :: sat
            Boolean checking if error occured                   :: err
            Print associated to error and unsatifiability       :: err_print
"""
# IP 11/03/19 how print documentation using Sphinx
# _____________________________________________________________________________
    generate_initial_state.__doc__
#    print generate_initial_state.__doc__
# _____________________________________________________________________________

    init_vars, err, err_print = algebraic_steady_state(params)
    if err > 0:
        if DEBUG > 2:
            print(err_print)
            assert(0)
        return [], False, err, err_print

    # if DEBUG > 2:
    #     print 'In generate_initial_state() params = ', params
    #     print 'len(params) = ', len(params)

    init_vars, sat, err, err_print = validate_and_improve_ss(init_vars,
                                                             params,
                                                             copy=False)
    if err > 0 or sat == False:
        if DEBUG > 2:
            print(err_print)
            assert(0)
        return [], sat, err, err_print

    sat, err, err_print = check_vars_constraints(init_vars, params)

    return init_vars, sat, err, err_print
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------


def algebraic_steady_state(params):
    """
    Compute the initial conditions using algebraic equqations on the parameters
    Model specific code in model.algebraic_steady_state(params)

    INPUT:
            Parameter set :: params

    OUTPUTS:
            Generated initial conditions    :: init_vars
            Boolean checking errors         :: err
            Print associated to errors      :: err_print
"""
# IP 11/03/19 how print documentation using Sphinx
# _____________________________________________________________________________
    algebraic_steady_state.__doc__
#    print algebraic_steady_state.__doc__
# _____________________________________________________________________________

    init_vars, err, err_print = model.algebraic_steady_state(params)

    return init_vars, err, err_print
# -----------------------------------------------------------------------------

# ----------------------------------------------------------------------------


def validate_and_improve_ss(init_vars, params, copy=True):
    """
    Check if an initial condition is a steady state.
    If not try 'MAX_ITERATION_SS_IMPROVEMENT' times to simulate and converge on
    the steady state

    INPUTS:
            Previous initial conditions to check    :: init_vars
            Parameter set                           :: params
            In-place or copy mode                   :: copy

    OUTPUTS:
            Steady state initial conditions                       :: init_vars
            Boolean checking if the initial conditions
            are steady states                                     :: sat
            Boolean checking errors                               :: err
            Print associated to errors or unsatisfiability         :: err_print
"""
# IP11/03/19 how print documentation using Sphinx
# _____________________________________________________________________________
    validate_and_improve_ss.__doc__
#    print validate_and_improve_ss.__doc__
# _____________________________________________________________________________

    if copy == True:
        init_vars = init_vars.copy()

    err = 0
    err_print = 'No error'
    sat = False

    # MAX_ITERATION_SS_IMPROVEMENT = 10

    nb_ite = 0
    while (nb_ite <= config.MAX_ITERATION_SS_IMPROVEMENT):
        nb_ite += 1
        simu_sols, err, err_print = simulate(model.SS_SIMU_ARGS,
                                             params,
                                             init_vars,
                                             config.SIMU_RTOL,
                                             config.SIMU_ATOL,
                                             config.Fe_replen_ss,
                                             False)
        if err > 0:
            if DEBUG > 0:
                print('Error in simulate() function -- numerical simulation')
            if DEBUG > 2:                                 # Anormal behaviour
                print(err_print)
                assert(0)
            return [], sat, err, err_print

        if DEBUG > 1:
            if nb_ite == 1:
                print('\nInitial state: ', simu_sols[0][0], '\n')

        # Get the time at which the system
        # reached a steady state (if not yet reached, simulate longer...)
        t_index_ss, err, err_print = check_stability(simu_sols,
                                                     params,
                                                     model.SS_SIMU_ARGS,
                                                     config.EPSILON)
        if err > 0:
            if DEBUG > 0:
                print('Error in check_stability().')
            if DEBUG > 2:  # Anormal behaviour
                print(err_print)
                assert(0)
            return [], sat, err, err_print

        var_sols = simu_sols[0]
        # Rq: t_index_ss = -1 take the last element
        init_vars[:] = var_sols[t_index_ss]

        if t_index_ss >= 0:
            sat = True
            if DEBUG > 1:
                print('Steady state found after ', nb_ite, 'iteration(s).')
                print('Initial state obtained by integrating the ODEs:')
                print('traces[', t_index_ss, '] =', var_sols[t_index_ss], '\n')

            # Normal behaviour output (steady state found)
            return init_vars, sat, err, err_print

        else:  # No steady state found Continue simulation (with Fe_ext = cste)

            if DEBUG > 1:
                print('    Steady state NOT FOUND, yet. Continuing...')
                print('_______________________________________________' + \
                    '______________________________')
            pass

    err = 0
    err_print = ' No Steady State found after {} iterations.'\
        .format(config.MAX_ITERATION_SS_IMPROVEMENT)

    return [], sat, err, err_print  # no steady state founds output only
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------


def simulate(simu_infos, params, init_vars, rtol, atol, inputs_list, perturbation = False):
    """
    Solve a set of differential equations using the solver defined in
    NUMERICAL_SOLVER.

    INPUTS:
            Information of the simulated system         :: simu_infos
                    -> pair (diff_eqs,t_range)
                       with diff_eqs a pointer on a function,
                       and t_range the list outputed time
            Parameter set                               :: params
            Initial conditions                          :: init_vars
            Relative tolerance                          :: rtol
            Absolute tolerance                          :: atol

    OUTPUTS:

            Simulation results          :: simu_sol
                    --> Tuple (trace_solutions, t_range, init_vars)
            Boolean checking errors     :: err
            Associated print            :: err_print
"""
# IP 11/03/19 how print documentation using Sphinx
# _____________________________________________________________________________
    simulate.__doc__
#    print simulate.__doc__
# _____________________________________________________________________________

    err = 0
    err_print = 'No error'

    if config.NUMERICAL_SOLVER == 'scipy-odeint-default':
        simu_sol, err, err_print = integrate_scipy_default(simu_infos,
                                                           params,
                                                           init_vars,
                                                           rtol,
                                                           atol,
                                                           inputs_list,
                                                           perturbation
                                                           )
        # some examples if one want to change the numerical solver
    elif config.NUMERICAL_SOLVER == 'scipy-solve_ivp-LSODA':
        simu_sol, err, err_print = integrate_scipy_solve_ivp(simu_infos,
                                                             params,
                                                             init_vars,
                                                             rtol,
                                                             atol,
                                                             inputs_list,
                                                             perturbation
                                                             )
        # TODO implementation see scipy version 1.0 parameteres may need to be
        # passed as global ? or see using lambda

        # err = 1
        # err_print = 'NUMERICAL_SOLVER {} not implemented.\n>>>> use default:
        # scipy-odeint-default '
        # err_print.format(NUMERICAL_SOLVER)
        # if DEBUG > 2:
        #    print(err_print)
        #    assert(0)
        # return [], err, err_print

    elif config.NUMERICAL_SOLVER == 'numba-odeint-default':
        err = 1
        err_print = 'NUMERICAL_SOLVER {} not implemented.\n>>>> use default:\
        scipy-odeint-default '
        err_print.format(config.NUMERICAL_SOLVER)
        if DEBUG > 2:
            print(err_print)
            assert(0)
        return [], err, err_print

        return [], err, err_print

    else:
        err = 1
        err_print = 'Numerical integration method {} not implemented.\n>>>> \
        use default: scipy-odeint-default '
        err_print.format(config.NUMERICAL_SOLVER)
        if DEBUG > 2:
            print(err_print)
            assert(0)
        return [], err, err_print
    return simu_sol, err, err_print
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------


def integrate_scipy_default(simu_infos, params, init_vars, rtol, atol, inputs_list, perturbation = False):
    """
    Solve some ODEs using odeint default solver from scipy

    INPUTS:
            Information of the simulated system         :: simu_infos
                    -> pair (diff_eqs,t_range)
                       with diff_eqs a pointer on a function,
                       and t_range the list outputed time
            Parameter set                               :: params
            Initial conditions                          :: init_vars
            Relative tolerance                          :: rtol
            Absolute tolerance                          :: atol

    OUTPUTS:
            Simulation results          :: simu_sol
                    --> Tuple (trace_solutions, t_range, init_vars)
            Boolean checking errors     :: err
            Associated print            :: err_print
"""
# IP 11/03/19 how print documentation using Sphinx
# _____________________________________________________________________________
    integrate_scipy_default.__doc__
#    print integrate_scipy_default.__doc__
# _____________________________________________________________________________

    err = 0
    err_print = 'No error'
    # 'pointer' on differential equation function (built at the import time)
    diff_eqs = simu_infos[0]
    t_range = simu_infos[1]

    try:
        var_sols = scipy.integrate.odeint(diff_eqs,
                                          init_vars,
                                          t_range,
                                          args=(params,inputs_list,perturbation),
                                          rtol=rtol,
                                          atol=atol)
        # var_sols = scipy.integrate.odeint(diff_eqs, init_vars,
        # t_range, rtol=rtol, atol=atol)

    except: # TODO bad exception handling
        err = 1
        err_print = 'In integrate_scipy_default():\n>>> Error in ' + \
            ' numerical simulation trying to solve \n ODE {} \n' + \
            '_______________________________________________' + \
            '_______________________________________________'

        err_print = err_print.format(diff_eqs)
        if DEBUG > 2:
            print(err_print)
            assert(0)
        return [], err, err_print

    simu_sol = (var_sols, t_range, init_vars)
    return simu_sol, err, err_print
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# TODO Still not working


def integrate_scipy_solve_ivp(simu_infos, params, init_vars, rtol, atol):

    err = 0
    err_print = 'No error'
    # 'pointer' on differential equation function (built at the import time)
    diff_eqs = simu_infos[0]
    t_range = simu_infos[1]

    try:
        # assert(False)
        sol = scipy.integrate.solve_ivp(lambda t,
                                        y: diff_eqs(y, t, params),
                                        t_range, init_vars,
                                        method='LSODA',
                                        rtol=rtol,
                                        atol=atol)
        # var_sols = scipy.integrate.odeint(diff_eqs, init_vars,
        # t_range, rtol=rtol, atol=atol)

    except:  # TODO bad exception handling
        err = 1
        err_print = 'In integrate_scipy_solve_ivp():\n>>> Error in numerical \
        simulation.\n trying to solve ODE {}'
        err_print = err_print.format(diff_eqs)
        if DEBUG > 2:
            print(err_print)
            assert(0)
        return [], err, err_print

    if sol.success == False:
        err = 1
        err_print = 'In integrate_scipy_solve_ivp():\n>>> Error in numerical \
        simulation.\n trying to solve ODE {}\n'
        err_print += sol.message
        return [], err, err_print

    simu_sol = (sol.y, sol.t, init_vars)
    return simu_sol, err, err_print
# -----------------------------------------------------------------------------

# ----------------------------------------------------------------------------


def check_stability(simu_sols, params, simu_infos, epsilon):
    """
    Check if a trajectory solution contain a steady state.

    INPUTS:
            Simulation results                      :: simu_sols
                    --> Tuple (trace_solutions, t_range, init_vars)
            Parameter set                           :: params
            Informations on the simulated system    :: simu_infos
                    -> pair (diff_eqs,t_range)
                       with diff_eqs a pointer on a function,
                       and t_range the list outputed time
            Maximum authorized variation            :: epsilon

    OUTPUTS:
            Time index of first entry in steady state   :: time_idx
            Boolean checking the error                  :: err
            Print associated to the error               :: err_print

"""
# IP 11/03/19 how print documentation using Sphinx
# _____________________________________________________________________________
    check_stability.__doc__
#    print check_stability.__doc__
# _____________________________________________________________________________

    var_sols, t_range, init_vars = simu_sols[0], simu_sols[1], simu_sols[2]
    diff_eqs = simu_infos[0]
    err = 0
    err_print = 'No error'
    nb_steady_it = 0

    for time_idx in range(len(var_sols) - 1):         # for each time point
        current_derivatives = None
        try:
            current_derivatives = diff_eqs(var_sols[time_idx],
                                           t_range[time_idx],
                                           params,
                                           config.Fe_replen_ss,
                                           False)
            # current_derivatives = diff_eqs(var_sols[time_idx],
            # t_range[time_idx])  # if params passed by global variables
        except:
            err = 1
            err_print = 'Error in computing derivatives: \n       >> At time \
            point idx {} deriv(t = {}) = {} and solutions = {}'
            err_print = err_print.format(time_idx,
                                         t_range[time_idx],
                                         current_derivatives,
                                         var_sols[time_idx])
            if DEBUG > 2:
                print(err_print)
                assert(0)
            return -1, err, err_print

        nvars = len(current_derivatives)

        nb_steady_vars = 0
        try:
            for ivar in range(nvars):  # for each var
                # normalized derivative < epsilon
                if abs(current_derivatives[ivar] / var_sols[
                        time_idx][ivar]) < epsilon:
                    nb_steady_vars += 1
        except:
            err = 1
            err_print = 'Error in check_stability(): \n      >>At time point \
            idx {} deriv(t = {}) = {} and solutions = {}'
            err_print = err_print.format(time_idx,
                                         t_range[time_idx],
                                         current_derivatives,
                                         var_sols[time_idx])
            if DEBUG > 2:
                print(err_print)
                assert(0)
            return -1, err, err_print

        if nb_steady_vars == nvars:
            # return  time_idx, err, err_print
            nb_steady_it += 1
        else:
            nb_steady_it = 0

        if nb_steady_it == config.NB_STEADY_POINTS:
            # print 'nvars = ', nvars
            # print 'current_derivatives[ivar] = ', current_derivatives[5]
            # print  'var_sols[time_idx][ivar] = ', var_sols[time_idx][5]
            # print 'abs(current_derivatives[ivar] / var_sols[time_idx][ivar])
            # = ', abs(current_derivatives[5] / var_sols[time_idx][5])
            return time_idx, err, err_print

    return -1, err, err_print
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------


def check_vars_constraints(init_vars, params):
    """
    Check absolute constraints on the initial conditions

    INPUTS:
            Initial conditions  :: init_vars
            Parameter set       :: params

    OUTPUTS
            Boolean checking the statisfiability of the
            constraints                                            :: sat
            Boolean checking if errors occured                     :: err
            Print associated to errors and unsatifiability         :: err_print
"""
# IP 11/03/19 how print documentation using Sphinx
# _____________________________________________________________________________
    check_vars_constraints.__doc__
#    print check_vars_constraints.__doc__
# _____________________________________________________________________________

    sat, err, err_print = model.check_initial_state_constraint(init_vars,
                                                               params)
    return sat, err, err_print
# ----------------------------------------------------------------------------

# ----------------------------------------------------------------------------
# TODO: The initial conditions differs from the one previously generated (one
#  additional variable)
#       -> either change the equations to remove this variable (Fe_ext) or add
#          a function which call the model.py module to update the initial
#          condition if necessary. This way we can avoid the model specific
#          code:  'init_vars+[config.Fe_GROWTH]'


def simulate_and_check(params, init_vars, inputs_list):
    """
    Simulates the model amd check absolute constraints. If there are multiple
    experiments, or input values, it performs mulitple simulations. ex: Iron
    re-apport in Iron Homeostasis model. After the simulation, it check results
    and compute additional informations using fonction provided in the model
    module.
    Currently there is no relative constraints check between simulations.
    The outputed solution is in a model dependant format >>> NOT ANYMORE

    INPUTS:
            Parameter set           :: params
            Initial conditions      :: init_vars
            List of possible Inputs :: inputs_list

    OUTPUTS:
            Solution of the multiple simulations                :: sol
                    -> Tuple (parameter set, relative solutions (ratios),
                             list of simulations results, initial conditions)
                                    -> Tuple: (trace, time range, initial
                                               conditions)

            Boolean checking the satisfiability of constraints  :: sat
            Boolean checking if errors occured                  :: err
            Print associated to errors or unsatifiability       :: err_print
"""
    all_simus = []
    sat = True
    for input_value in inputs_list:
        simu_sols, err, err_print = simulate(model.MAIN_SIMU_ARGS,
                                             params,
                                             init_vars,
                                             config.SIMU_RTOL,
                                             config.SIMU_ATOL,
                                             input_value,
                                             True
                                             )
#        # TODO model dependant init_vars+[config.Fe_GROWTH]
#        simu_sols, err, err_print = simulate(model.MAIN_SIMU_ARGS,
#                                             params+[input_value],
#                                             init_vars+[config.Fe_GROWTH],
#                                             config.SIMU_RTOL, config.SIMU_ATOL)
        # assert(False)

        if err > 0:
            return [], False, err, err_print

        sat, err, err_print = check_absolute_constraints(simu_sols,
                                                         input_value)
        if err > 0 and DEBUG > 0:
            print('params = ', params)
            print('input_value = ', input_value)
            print('init_vars = ', init_vars+[config.Fe_GROWTH])
            print('Last point = ', simu_sols[0][model.INDEX_END])
            if DEBUG > 1:
                print('\n Press Return to continue')
                input()
        if sat == False:
            return [], sat, err, err_print

        all_simus.append(simu_sols)

    relative_sols, err, err_print = compute_relative_values(all_simus)
    if err > 0:
        return [], False, err, err_print

    sat, err, err_print = check_relative_constraints(relative_sols)
    if sat == False:
        return [], sat, err, err_print

    sol = [params, relative_sols, all_simus, init_vars]

    return sol, sat, err, err_print
# ----------------------------------------------------------------------------

# ----------------------------------------------------------------------------


def check_absolute_constraints(simulation, current_input):
    """
    Check absolute constraints on the solution of the simulations.
    Actual model depedant code in model.check_absolute_constraints()

    INPUTS:
            Simulations results                                :: simulation
                -> Tuple: (trace, time range, initial conditions)
            Input used for this simulation                     :: current_input

    OUTPUTS
            Boolean checking the satifiability of the constraints  :: sat
            Boolean checking if errors occured                     :: err
            Print associated to errors or unsatifiablity           :: err_print
"""
# IP 11/03/19 how print documentation using Sphinx
# _____________________________________________________________________________
    check_absolute_constraints.__doc__
#    print check_absolute_constraints.__doc__
# _____________________________________________________________________________

    sat, err, err_print = model.check_absolute_constraints(simulation,
                                                           current_input)

    return sat, err, err_print
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------


def compute_relative_values(simulations_list):
    """
    Compute relative values between the results of multiple simulations
    Actual model specific code in model.compute_relative_values()

    INPUT:
            List of simulation results  :: simulations_list

    OUTPUTS:
            Compute relative solutions              :: relative sols
                    -> List of tuples [ (ratio1_input1, ratio2_input1, ...),
                                        (ratio1_input2, ratio2_input2, ...),
                                        ... ]
            Boolean checking if errors occured      :: err
            Print associated to errors              :: err_print

"""
# IP 11/03/19 how print documentation using Sphinx
# _____________________________________________________________________________
    compute_relative_values.__doc__
#    print compute_relative_values.__doc__
# _____________________________________________________________________________

    relative_sols, err, err_print =\
        model.compute_relative_values(simulations_list)

    return relative_sols, err, err_print
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------


def check_relative_constraints(relative_result_list):
    """
        Check constraint on relative values.
        -> Do nothing
"""
# IP 11/03/19 how print documentation using Sphinx
# _____________________________________________________________________________
    check_relative_constraints.__doc__
#    print check_relative_constraints.__doc__
# _____________________________________________________________________________

    sat, err, err_print =\
        model.check_relative_constraints(relative_result_list)

    return sat, err, err_print
# -----------------------------------------------------------------------------


# -----------------------------------------------------------------------------
# TODO ERROR CODE IS NOT FINISHED


def print_end_infos(list_sols, file_name):
    """
     Print results in file on console.
     Default console output: number of sols and enclosing box
     If specified in model print additional infos

    INPUT:
            list of solutions obtained after all the simulations   :: list_sols

    OUTPUT:
            Boolean checking errors     :: err
            Associated print            :: err_print
"""
# IP 11/03/19 how print documentation using Sphinx
# _____________________________________________________________________________
    print_end_infos.__doc__
#    print print_end_infos.__doc__
# _____________________________________________________________________________

    err = 0
    err_print = 'ERROR NOT HANDLED'
    longest_string = 0
    for item in config.PARAM_NAMES:
        if len(item) > longest_string:
            longest_string = len(item)
# IP 08/03/19
# I changed the way to print the param bounds new in order to avoid 3 loops
# I directy zipped param bounds (from list_sols) and thier name from model and
# I sorted the values of param bounds in row using .value.sort() function

    if len(list_sols) >=1:
        print('\n=========== Simulate Found ', len(list_sols), 'solutions ' + \
            '============\n')
#        print 'number of parameters: ', N_par
        # IP 12/03/19
        # FORCED_SAMPLING_SCALE == 'linear' return a tupla list_sols which
        # constists in a n=len(list_sols) sets of 4 objects:
        # list_sols=[[params_i[p1,...p25, input=0]
        #            [relative_sols 9 elements]
        #            [all_simus is empty elements for memory reason:
        #            [init_vars - 7 elements]]
        list_par = [col[0] for col in list_sols]
        list_var = [col[3] for col in list_sols]
        list_rel_sol_i = [col[1] for col in list_sols]
        list_rel_sol = []
        for i in range(0,len(list_sols)):
            for j in range(0, len(list_sols[0][1])):
                list_rel_sol.append(list_rel_sol_i[i][j])
        list_FeTot = [col[2] for col in list_sols]
        list_parT = list(zip(*list_par))
        list_varT = list(zip(*list_var))
        list_rel_solT = list(zip(*list_rel_sol))
        list_FeTotT = list(zip(*list_FeTot))
        list_tot = list_parT + list_varT + list_rel_solT + list_FeTotT
        ratio_name = ['ratio_Fe-sto_250', 'ratio_TfR1_250', 'ratio_Ft_250',
                      'ratio_Fe-sto_25', 'ratio_TfR1_25', 'ratio_Ft_25',
                      'ratio_Fe-sto_10', 'ratio_TfR1_10', 'ratio_Ft_10']
        if len(list_rel_solT) == 0:
            ratio_name = []
        FeTot_name = ['Fe_Tot_250', 'Fe_Tot_25', 'Fe_Tot_10', 'Fe_Tot_0']
        if len(list_FeTot[0]) == 1:
            FeTot_name = ['Fe_Tot_250']
        indice = config.PARAM_NAMES + config.VAR_NAMES + ratio_name + FeTot_name
        name_col = ['set_' + str(jj) for jj in range(1, len(list_sols) + 1)]
        list_parT = pd.DataFrame(
                list_parT, index=config.PARAM_NAMES, columns=name_col)
        list_varT = pd.DataFrame(
                list_varT, index=config.VAR_NAMES, columns=name_col)
        list_rel_solT = pd.DataFrame(
                list_rel_solT, index=ratio_name, columns=name_col)
        list_FeTotT = pd.DataFrame(
                list_FeTotT, index=FeTot_name, columns=name_col)
        list_tot = pd.DataFrame(
                list_tot, index=indice, columns=name_col)
        # IP 03 04 19 write results of parameters set into file txt
        # save into file before sort
        results_file = os.path.join(run_path, rfile)
        with open(results_file, 'a') as fileR:
            fileR.write(list_tot.to_string() + '\n')
        # Free large summary table ASAP; keep only sorted param bounds for print.
        del list_tot, list_varT, list_rel_solT, list_FeTotT
        list_parT.values.sort()
        if len(list_sols) > 1:
            tab_par_b_sols = list(zip(
                    list_parT['set_1'],
                    list_parT['set_' + str(len(list_sols))],
                    config.PARAM_NAMES))
            tab_par_bounds_sols = pd.DataFrame(
                    tab_par_b_sols, np.arange(1, config.NB_DIM_TOT+1), columns=[
                            'lower bounds', 'upper bounds', 'name'])
            print('============ Parameters bounds table ===============')
            print(tab_par_bounds_sols)
            print('====================================================')
        elif len(list_sols) == 1:
            # IP 16/07/19 write result of only one set into file txt
            list_par = [list_sols[0][0]]
            list_parT = list(itertools.chain.from_iterable(list_par))
            list_parT = pd.DataFrame(
                    list_parT, index=config.PARAM_NAMES, columns=['set_1'])
            print('============ Parameters bounds table ===============')
            print(list_parT)
            print('====================================================')
    else:
        print('\n =========== Simulate NO solutions found ============')

    return err, err_print  # TODO err handling
# -----------------------------------------------------------------------------


# -----------------------------------------------------------------------------
# TODO re-read and complete DEBUG and error


def produce_over_box(ifile_name, size_incr):
    """
    Description

    Inputs: String of the input file name   :: ifile_name
            String of the output file name  :: ofile_name
            Size increase (percentage)      :: size_incr

    Ouput:
"""
# IP 11/03/19 how print documentation using Sphinx
# _____________________________________________________________________________
    produce_over_box.__doc__
#    print produce_over_box.__doc__
# _____________________________________________________________________________

    # Read the parameter bounds in the file
    # (each line -> bounds for 1 parameter)
    if DEBUG > 1:
        print('reading file', ifile_name)
    # Note: the values are indexed according to their order in the input file

    param_bounds, nb_params, sampled_names = read_file.read_file(
            ifile_name, DEBUG, config.map_PARAM_NAMES)

    ofile_name = 'inputs/overbox-'+str(size_incr)+'percent-'+ifile_name[9:]
    try:
        ofile = open(ofile_name, 'w')
    except IOError:
        func_name = inspect.currentframe().f_code.co_name
        print('\n*** ERROR in %s: cannot open file %s\nabort...\
            ' % (func_name, ifile_name), file=sys.stderr)
        assert(False)

    try:
        longest = 0
        for item in sampled_names:
            if len(item) > longest:
                longest = len(item)
        fmt_print_param = '{vname: <{width}} :  {min}    {max} \n'

        for ipar in range(nb_params):
            print('param_bounds[', ipar, '] = ', param_bounds[ipar])
            ext_min = param_bounds[ipar][0] - (
                    size_incr/100.0)*param_bounds[ipar][0]
            ext_max = param_bounds[ipar][1] + (
                    size_incr/100.0)*param_bounds[ipar][1]
            ofile.write(fmt_print_param.format(
                    vname=sampled_names[ipar],
                    width=longest,
                    min=ext_min,
                    max=ext_max))
    finally:
        ofile.close()
    return
# -----------------------------------------------------------------------------


# -----------------------------------------------------------------------------
def _empty_trajectory_store():
    """Nested lists: [n_vars][n_inputs] -> list of 1-D trajectory arrays."""
    n_var = len(config.VAR_NAMES)
    n_in = len(config.MULTIPLE_INPUTS)
    return [[[] for _ in range(n_in)] for _ in range(n_var)]


def _append_trajectories(store, sol):
    """Copy per-variable columns out of sol so the full ODE buffer can be GC'd."""
    for ivar in range(len(config.VAR_NAMES)):
        for iin in range(len(config.MULTIPLE_INPUTS)):
            store[ivar][iin].append(np.asarray(sol[2][iin][0][:, ivar]).copy())


def calculate_sol_y_t(list_sols):
    solutions = _empty_trajectory_store()
    N_sol = len(list_sols)
    if N_sol < 1:
        return solutions

    N_set = max(1, int(config.MAX_TRAJECTORIES_DUMP))
    pace = 1
    if N_sol >= N_set:
        pace = int(N_sol / N_set)
    for isol in range(0, N_sol, pace):
        sol, sat, err, err_print = simulate_and_check(
                list_sols[isol][0], list_sols[isol][3], config.MULTIPLE_INPUTS)
        _append_trajectories(solutions, sol)
        del sol
    return solutions
# -----------------------------------------------------------------------------


# -----------------------------------------------------------------------------
def calculate_sol_filtered_y_t(list_sols, saved_index_sol):
    solutions_filtered = _empty_trajectory_store()
    N_sol = len(saved_index_sol)
    if N_sol < 1:
        print('The solutions files are not written.')
        return solutions_filtered

    N_set = max(1, int(config.MAX_TRAJECTORIES_DUMP))
    if N_sol > N_set:
        stride = max(1, int(N_sol / N_set))
        for index, isol in enumerate(saved_index_sol, start=1):
            if index % stride != 0:
                continue
            sol, sat, err, err_print = simulate_and_check(
                    list_sols[isol][0], list_sols[isol][3],
                    config.MULTIPLE_INPUTS)
            _append_trajectories(solutions_filtered, sol)
            del sol
    else:
        for isol in saved_index_sol:
            sol, sat, err, err_print = simulate_and_check(
                    list_sols[isol][0], list_sols[isol][3],
                    config.MULTIPLE_INPUTS)
            _append_trajectories(solutions_filtered, sol)
            del sol
    return solutions_filtered
# -----------------------------------------------------------------------------
def dump_solution_into_files(solutions, path_sol):
    # IP need to pass calculate_sol_y_t or calculate_filter_sol_y_t
    N_sol = len(solutions[0][0])
    if N_sol < 1:
        print('The solutions files are not written.')
        return
    fe_replen = ['_250_', '_25_', '_10_', '_0_']

    for i in range(0, len(config.VAR_NAMES)):
        for kk in range(0, len(config.MULTIPLE_INPUTS)):
            results_file = os.path.join(
                    path_sol, 'results_y' + config.VAR_NAMES[i] +
                    fe_replen[kk]+'.txt')
            data_print = []
            names = [config.VAR_NAMES[i] + '_' + fe_replen[kk] +
                     str(m) for m in range(1, N_sol + 1)]
            col = pd.Index(names, name='time')
            for jj in range(0, N_sol):
                data_print.append(np.asarray(solutions[i][kk][jj]))
            data_printT = list(zip(*data_print))
            with open(results_file, 'a') as fileR:
                fileR.write(pd.DataFrame(data_printT,
                                         index=g_time,
                                         columns=col).to_string() +
                                         '\n')
    return
# ----------------------------------------------------------------------------


# -----------------------------------------------------------------------------
def plot_results(solutions, plot_path):
    N_sol = len(solutions[0][0])
    x = g_time / 3600.0
    fe_replen = ['_250_', '_25_', '_10_', '_0_']
    if N_sol >= 1:
        for i in range(0, len(config.VAR_NAMES)):
            for kk in range(0, len(config.MULTIPLE_INPUTS)):
                fig = plt.figure()
                plt.xlabel('time (h)')
                plt.title(config.VAR_NAMES[i] +
                          ' Concentration (M) at Fe replen =' +
                          str(config.MULTIPLE_INPUTS[kk]) + ' M')
                plt.ylabel(config.VAR_NAMES[i] +
                           ' Concentration (M)')
                pace = 1
                N_curves = max(1, int(config.MAX_PLOT_CURVES))
                if N_sol >= N_curves:
                    pace = int(N_sol / N_curves)
                for jj in range(0, N_sol, pace):
                    y = np.asarray(solutions[i][kk][jj])
                    plt.plot(x, y)
                plt.savefig(plot_path + '/' + 'Plot_y' + config.VAR_NAMES[i] +
                            fe_replen[kk] + '.png', dpi=200)
                plt.close(fig)
    else:
        print('The plot file results_y' + config.VAR_NAMES[0] +
              '.txt is not written.')
    return
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
def plot_correlation_param(list_sols):
    list_par = [col[0] for col in list_sols]
    data_param = np.array(list_par)
    data_paramT = list(zip(*data_param))
    for ii in range(0, len(data_paramT) - 1):
        for jj in range(ii + 1, len(data_paramT)):
            fig = plt.figure()
            plt.scatter(data_paramT[ii], data_paramT[jj], marker='.')
            plt.xlim([min(data_paramT[ii]) * 0.1, max(data_paramT[ii]) * 10])
            plt.ylim([min(data_paramT[jj]) * 0.1, max(data_paramT[jj]) * 10])
            plt.xlabel(config.PARAM_NAMES[ii])
            plt.ylabel(config.PARAM_NAMES[jj])
            plt.yscale('log')
            plt.xscale('log')
            plt.savefig(path_plot_correlation_param1 +
                        '/' + str(ii) + '-' + str(jj) + 'Plot_' +
                        config.PARAM_NAMES[ii] + 'VS' +
                        config.PARAM_NAMES[jj] + '.png', dpi=200)
            plt.close(fig)
    return
# -----------------------------------------------------------------------------


# -----------------------------------------------------------------------------
def plot_correlation_param_filter(list_sols, saved_index_sol):
    list_par = [col[0] for col in list_sols]
    data_param = np.array(list_par)
    data_param_filter = data_param[saved_index_sol]
    data_param_filterT = list(zip(*data_param_filter))
    for ii in range(0, len(data_param_filterT) - 1):
        for jj in range(ii + 1, len(data_param_filterT)):
            fig = plt.figure()
            plt.scatter(data_param_filterT[ii],
                        data_param_filterT[jj], marker='.')
            plt.xlim([min(data_param_filterT[ii]) -
                      0.15 * min(data_param_filterT[ii]),
                      max(data_param_filterT[ii]) +
                      0.05 * max(data_param_filterT[ii])])
            plt.ylim([min(data_param_filterT[jj]) -
                      0.15 * min(data_param_filterT[jj]),
                      max(data_param_filterT[jj]) +
                      0.05 * max(data_param_filterT[jj])])
            plt.xlabel(config.PARAM_NAMES[ii])
            plt.ylabel(config.PARAM_NAMES[jj])
            plt.yscale('log')
            plt.xscale('log')
            plt.savefig(path_plot_correlation_param + '/' + str(ii) + '-' +
                        str(jj) + 'Plot_' +
                        config.PARAM_NAMES[ii] + 'VS' +
                        config.PARAM_NAMES[jj] + '.png', dpi=200)
            plt.close(fig)
    return
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
def compute_Fe_tot(all_simus):
    Fetot = [0]*config.NB_INPUTS
    for i in range(config.NB_INPUTS):
        traces = all_simus[i][0]
        # Currently INDEX_END correspond to t = 53h as t_cutoff ~ 5h
        Fe = traces[model.INDEX_END][0]
        # Ft:Fe[end of re-entry exp]
        FtFe = traces[model.INDEX_END][2]
        Festo = traces[model.INDEX_END][3]
        # (Fe + Ft:Fe + Fe-sto)[end of re-entry exp]
        Fetot[i] = Fe + FtFe + Festo
    return Fetot
# -----------------------------------------------------------------------------


def _stable_err_print(err_print):
    """Drop volatile addresses from error strings (function repr changes)."""
    text = str(err_print)
    # e.g. <function diff_eqs at 0x000001FE...>
    import re
    return re.sub(r' at 0x[0-9A-Fa-f]+', '', text)


def evaluate_parameter_set(sampled_params, inputs_list=None):
    """
    Run one sampled parameter vector through the HIF pipeline and return a
    JSON-serializable fingerprint used by regression tests.

    sampled_params: sequence of length NB_SAMPLED_DIM (or NB_DIM_TOT).
    """
    if inputs_list is None:
        inputs_list = config.MULTIPLE_INPUTS

    fingerprint = {
        'sat': False,
        'err': 0,
        'stage': 'start',
        'err_print': '',
        'algebraic_init_vars': [],
        'init_vars': [],
        'end_states': [],
        'relative_sols': [],
        'Fe_tot': [],
        'params': [],
    }

    params = list(sampled_params)
    err, err_print = check_negativity(params)
    if err > 0:
        fingerprint.update(
            stage='negativity', err=err,
            err_print=_stable_err_print(err_print))
        return fingerprint

    params, err, err_print = convert_to_recordclass(params)
    if err > 0:
        fingerprint.update(
            stage='convert', err=err,
            err_print=_stable_err_print(err_print))
        return fingerprint

    params, err, err_print = add_sub_sample(
        params, config.DEFAULT_ADD_PARAMS[0], copy=False)
    if err:
        fingerprint.update(
            stage='add_sub_sample', err=err,
            err_print=_stable_err_print(err_print))
        return fingerprint

    params, err, err_print = add_params(params, copy=False)
    if err:
        fingerprint.update(
            stage='add_params', err=err,
            err_print=_stable_err_print(err_print))
        return fingerprint

    params, sat, err, err_print = check_params(params, copy=False)
    if err or not sat:
        fingerprint.update(
            stage='check_params', sat=bool(sat), err=err,
            err_print=_stable_err_print(err_print),
            params=[float(params[i]) for i in range(config.NB_DIM_TOT)])
        return fingerprint

    fingerprint['params'] = [float(params[i]) for i in range(config.NB_DIM_TOT)]

    # Freeze algebraic steady-state even if numerical SS improvement fails.
    alg_init, alg_err, alg_err_print = algebraic_steady_state(params)
    if alg_err == 0 and alg_init:
        fingerprint['algebraic_init_vars'] = [float(x) for x in alg_init]

    init_vars, sat, err, err_print = generate_initial_state(params)
    if err or not sat:
        fingerprint.update(
            stage='initial_state', sat=bool(sat), err=err,
            err_print=_stable_err_print(err_print))
        return fingerprint

    sol, sat, err, err_print = simulate_and_check(params, init_vars, inputs_list)
    if err or not sat or not sol:
        fingerprint.update(
            stage='simulate_and_check', sat=bool(sat), err=err,
            err_print=_stable_err_print(err_print),
            init_vars=[float(x) for x in init_vars])
        return fingerprint

    end_states = []
    for simu in sol[2]:
        traces = simu[0]
        end_states.append([float(v) for v in traces[model.INDEX_END]])

    relative = []
    for item in sol[1]:
        if isinstance(item, (list, tuple)):
            relative.append([float(x) for x in item])
        else:
            relative.append(float(item))

    fingerprint.update(
        sat=True,
        err=0,
        stage='ok',
        err_print='',
        init_vars=[float(x) for x in init_vars],
        end_states=end_states,
        relative_sols=relative,
        Fe_tot=[float(x) for x in compute_Fe_tot(sol[2])],
        params=[float(params[i]) for i in range(config.NB_DIM_TOT)],
    )
    return fingerprint
# -----------------------------------------------------------------------------


# -----------------------------------------------------------------------------
def main():
    """
    Main fonction.
    call with:
    python [name_script.py]
           [var-bound-file-name]
           [par-bounds-file-name]
           [results-file-name]
           [MODE](optional)
           [NB-points]

    In python command line it can be launch using for example:
import simulate
ifolder = inputs/
ofolder = run_path
vfile = ifolder+ 'bounds-vars-initiales.list'
pfile = ifolder+ 'bounds-params-initiales.list'
rfile = ofolder+'temp-test-results.list'
npts = 10
simulate.apply_method( pfile, vfile, rfile, npts, MODE='Box-Sampling')

    Takes 4 madatory arguments:
       - The NAME of the file containing the bounds on the variables
       - The NAME of the file containing the bounds on the parameters OR
         a list of parameter points
       - The NAME of the result file (WILL be in outputs/)
                --> do not put the 'outputs/' in the file name as it is
                    hardcoded in the main.
       - The number of points in the simulation
    And 1 optional argument, the mode of the method ('Box-Sampling' by default)
    among:
        - 'Box-Sampling'
        - 'Plot'
        - 'Additional-Sampling'
        - 'Re-Sampling'
        - 'PCA-guided-sampling'
    Confers to function apply_method() for more details on the mode
"""
# IP 11/03/19 how print documentation using Sphinx
# _____________________________________________________________________________
    main.__doc__
#    print main.__doc__
# _____________________________________________________________________________
# IP 13_03_19 sys.argv define the option to run Simulate Script
# By default the 0 element is the name of the script

    start = time.time()
    model.DEBUG = DEBUG

    print('============================================== SIMULATE SCRIPT ' + \
          '==============================================')
    print('The variables bounds are read in ', vfile)
    print('The parameters bounds are read in ', pfile)
    print('The results will be write in ', rfile)

    # Inform the user of an inconsistency in the input as soon as possible
    if config.NB_SAMPLED_DIM + len(config.DEFAULT_ADD_PARAMS[0]) + \
            len(config.PARAM_ALGEB_DEFS) != config.NB_DIM_TOT:
        print('\nError in main(): the values of NB_SAMPLED_DIM,\
        DEFAULT_ADD_PARAMS, config.PARAM_ALGEB_DEFS and NB_DIM_TOT are\
        inconsistent')
        sys.exit(1)

    if (len(args_list) >= 6):
        mode = args_list[4]  # [option]
        npts = int(args_list[5])  # [NB-points]
        N_sols, list_sols = apply_method(pfile, vfile, rfile, npts, mode)
        # IP 14 03 19
        # after run sinulate, we call filter which use as input file the output
        # file rfile and will be write in run_path
        # N_sol = list_sols, but is defined as global
        # I should import also the MODE== Box-Sampling
        # from apply_method import MODE
        end = time.time()
        print('simulate time execution: ' + str(end - start) + '\n')

        if N_sols > 0:
            print('Simulate is writing the solutions into files.' +
                  ' It will take time...')
            solutions = calculate_sol_y_t(list_sols)
            dump_solution_into_files(solutions, path_solutions)
            print('Simulate is plotting and saving the solutions into files.' +
                  ' It will take time...')
            plot_results(solutions, plot_solutions)
            # Free trajectory buffers before correlation / filter work.
            del solutions
            plot_correlation_param(list_sols)
            filename = rfile
            print('Simulate filter the solutions')
            save_index_sol, list_solutions = Filter.Filter(filename, run_path)
            index = save_index_sol
            solutions_filtered = calculate_sol_filtered_y_t(list_sols, index)
            log.log_file(
                    run_path,
                    script,
                    vfile,
                    pfile,
                    rfile,
                    mode,
                    npts,
                    config.FORCED_SAMPLING_SCALE,
                    config.RESAMPLING_SIZE,
                    config.skip,
                    config.magnitude,
                    config.increment,
                    config.MAX_ITERATION_SS_IMPROVEMENT,
                    config.SIMU_RTOL,
                    config.SIMU_ATOL,
                    config.EPSILON,
                    config.NB_STEADY_POINTS,
                    config.NUMERICAL_SOLVER,
                    config.Name_Model,
                    config.VAR_NAMES,
                    config.PARAM_NAMES,
                    config.PARAM_ALGEB_DEFS,
                    config.NB_DIM_TOT,
                    config.NB_SAMPLED_DIM,
                    config.DEFAULT_ADD_PARAMS,
                    config.NB_VARS,
                    config.MULTIPLE_INPUTS,
                    config.Fe_GROWTH,
                    config.VARS_BOUNDS,
                    config.PARAMS_CONSTRAINT_BOX,
                    config.Fe_By_Ft_CONSTRAINT,
                    config.MAX_Fe_CONCENTRATION,
                    config.CSTRT_Fetot_HIGH_END,
                    config.CSTRT_Fetot_LOW_END,
                    config.MIN_FREE_IRON,
                    config.NB_RATIOS,
                    config.START,
                    config.INC_m,
                    config.END_SS_m,
                    config.END_MAIN_m,
                    config.constraints_on_Fe_sto,
                    config.constraints_on_TfR1,
                    config.constraints_on_Ft,
                    config.Constraints_on_Ft0,
                    config.Constraints_on_Fetot,
                    config.upperbound_on_C,
                    list_sols,
                    index,
                    list_solutions)
            if len(save_index_sol) > 0:
                print('Simulate is writing the solutions into files.' +
                      ' It will takes time...')
                dump_solution_into_files(
                        solutions_filtered, path_filter_solutions)
                print('Simulate is ploting and saving the solutions into' +
                      ' files. It will takes time...')
                plot_results(solutions_filtered, plot_filter_solutions)
                plot_correlation_param_filter(list_sols, index)
            del solutions_filtered

        else:
            sys.exit(1)
    else:
        npts = int(args_list[4])  # [NB-points]
        apply_method(pfile, vfile, rfile, npts)
    return

# -----------------------------------------------------------------------------


if __name__ == '__main__':
    setup_run_environment()
    main()
