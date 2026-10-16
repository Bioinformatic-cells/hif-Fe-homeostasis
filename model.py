# -*- coding: utf-8 -*-
"""
Created on Mon Mar 30 14:21:00 2020

@author: irina
"""

# ! /usr/bin/env python
# -*- coding: utf-8 -*-

"""
Module containing all the functions and constant specific to a model.
"""

import numpy as np
from math import isnan, isinf, sqrt
import scipy.linalg as la
import pandas as pd
import os
# _________constant variable of model.py _____________________________________
from config import VAR_NAMES, PARAM_NAMES,TIME_CUTOFF
from config import NB_DIM_TOT, NB_SAMPLED_DIM
from config import MULTIPLE_INPUTS, Fe_GROWTH, O2_INPUT
from config import PHD_tot, HIF_PHD, O2_hypoxia, HIF_INIT_CONSTRAINT
from config import VARS_BOUNDS, PARAMS_CONSTRAINT_BOX, FT_INIT_CONSTRAINT
from config import Fe_By_Ft_CONSTRAINT, MAX_Fe_CONCENTRATION
from config import CSTRT_Fetot_HIGH_END, CSTRT_Fetot_LOW_END, MIN_FREE_IRON
from config import NB_RATIOS, DEBUG_M as DEBUG, MAX_HIF_CONCENTRATION
from config import INC, END_SS, END_MAIN, START,MAX_FtFe_CONCENTRATION

import warnings
warnings.filterwarnings("error")
iteration=0

# Stiffness log handle. Kept optional so the model can be imported without
# running the full simulate_* CLI (needed for regression tests).
class _NullStiffFile(object):
    def write(self, *args, **kwargs):
        return None

    def flush(self):
        return None

    def close(self):
        return None

fileStiff = _NullStiffFile()
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
"""List of the constants used for this model.
The constants stated as [MANDATORY] are used in the method module
(simulate.py), and should be defined.
ADD_SAMPLE_FILE         :: [MANDATORY] File containing bounds on the additional
                                       parameters
NB_ADD_SAMPLES          :: [MANDATORY] Number of samples when sampling the
                                       additional parameters
RESAMPLING_SIZE         :: [MANDATORY] % size to define the box of the
                                       re-sampling method
PCA_REDUCTION           :: [MANDATORY] Flag defining if there is dimension
                                      reduction or not in PCA method
                                      Accepted values are:  'reduced' or 'full'
PARAM_NAMES             :: [MANDATORY] Names and order of the parameters (used
                                       for printing only)
PARAM_ALGEB_DEFS        :: List of parameters which are algebraically defined
                           in add_params()
n                       :: Sigmoids degree in differential equations
NB_DIM_TOT              :: [MANDATORY] Total dimension of the parameters space
                                       (sampled and not sampled parameters)
                           TODO -> Could be obtained by PARAM_NAMES length
NB_SAMPLED_DIM          :: [MANDATORY] Dimension of the sampled parameter space
DEFAULT_ADD_PARAMS_NAMES:: List of the parameters names for the parameters with
                           values fixed by default
DEFAULT_ADD_PARAMS      :: [MANDATORY] list of default values for fixed
                                       parameters, '[[]]' if there is none.
VAR_NAMES               :: [MANDATORY] List of the variables names (whose
                                       initial state is generated)
NB_VARS                 :: [MANDATORY] Number of variables (whose initial state
                                       is generated)
                           TODO -> Could be obtained using VAR_NAMES length
MULTIPLE_INPUTS         :: [MANDATORY] Possible values of the inputs
                                       (at least 1)
Fe_GROWTH               :: Value of the external iron concentration in the
                           Growing case
TIME_CUTOFF             :: Effective cutoff time of the external iron
VARS_BOUNDS             :: Constraint bounds on the initial conditions
PARAMS_CONSTRAINT_BOX   :: Constraint bounds on the sampled parameters
FT_INIT_CONSTRAINT      :: (DUPLICATED OF ABOVE) Constraint on the initial
                           value of Ft
Fe_By_Ft_CONSTRAINT     :: Maximum iron load by Ferritin
MAX_Fe_CONCENTRATION    :: Maximum iron concentration during the whole
                           simulation
CSTRT_Fetot_HIGH_END    :: Upper bounds on the [Total Iron - residual Iron] at
                           the end ---> Fe-FIC + Fe-sto + Ft:Fe
CSTRT_Fetot_LOW_END     :: Lower bounds on the [Total Iron - residual Iron] at
                           the end
MIN_FREE_IRON           :: Lower bounds on the Iron at the end of iron
                           deprivation
NB_RATIOS               :: Number of ratios for each possible replenishment
                           value (ie for each experiment)
DEBUG                   :: Debug mode:
                                       -  0 -> No trace
                                       -  1 -> error output
                                       -  2 -> more output
                                       -  3 -> even more and possibility of
                                               stopping the simulation
                           TODO -> Should be set by the simulate.py module
"""
# -----------------------------------------------------------------------------


# -----------------------------------------------------------------------------
# TODO: Note that the 2nd argument 't' is not used in the function.
"""Define differential system corresponding to Fe_ext constant: Steady-state
 case

    INPUTS:
            variable values :: y
            time            :: t
            parameters set   :: p
                C                  = p[0]
                k_Aco1_IRP1        = p[1]
                k_Fe_cons          = p[2]
                k_Fe_export        = p[3]
                k_Fe_input         = p[4]
                k_Ft_deg           = p[5]
                k_Ft_prod          = p[6]
                k_IRP1_Aco1        = p[7]
                k_IRP1_deg         = p[8]
                k_IRP_TfR1         = p[9]
                k_TfR1_deg         = p[10]
                k_TfR1_prod        = p[11]
                k_aso              = p[12]
                k_cat              = p[13]
                k_dis              = p[14]
                theta_Aco1_IRP1    = p[15]
                theta_IRP_Ft       = p[16]
                theta_IRP_TfR1     = p[17]
                k_FPN1a_deg        = p[18]
                k_FPN1a_prod       = p[19]
                k_IRP_FPN1a        = p[20]
                theta_IRP_FPN1a    = p[21]
                k_hif_prod         = p[22]
                k_hif_deg          = p[23]
                k_pdh_cat          = p[24]
                k_phd_dis          = p[25]
                k_phd_M            = p[26]
                theta_HIF_TfR1     = p[27]
                k_HIF_TfR1         = p[28]
                theta_IRP_HIF      = p[29]
                k_IRP_HIF          = p[30]

                k_IRP_Ft           = p[31]
                theta_IRP1_Aco1    = p[32]
                k_IRP_prod         = p[33]
    OUTPUT:
            set of derivatives :: ydot
"""
n = 2
min_ydot = np.zeros(7, dtype=float)
max_ydot = np.zeros(7, dtype=float)

ratios_max_min = []

def diff_eqs(y, t, p, inputs_list, perturbation):
    # IP 10/11/19 I introduced perturbation argument.
    # perturbation = False in the s.s., and True in the experiments.
    # IP 10/12/19 I removed input. I pass the multiple inputs as args of
    # diff_eqs()
    # Global variables used: t_cutoff, and parameters  (not concentrations)
    ydot = np.zeros((len(y)), dtype=float)    # create ODE vector

    Fe, Ft, FtFe, Fe_sto, IRP1, TfR1, HIF  = y[0], y[1], y[2], y[3], y[4],\
        y[5], y[6]  # Without Fe_ext as variable
    speed_kasoFeFt = p[12]*Fe*Ft
    speed_kFtdegFesto = p[5]*Fe_sto
    speed_kdispkcat_FtFe = (p[14] + p[13])*FtFe
    
    global min_ydot 
    global min_ydot
    global ratios_max_min

    # Fe_ext
    #Fe_ext = Fe_GROWTH
    
    # IP 4/10/19 needs to be changed for the new simulation
    sigm_time = 1. / ((t / (TIME_CUTOFF*3600))**20 + 1.)
    if t < 1e-30:
        sigm2_time = 0.
    else:
        sigm2_time = 1. / ((t / ((TIME_CUTOFF + 70)*3600))**20 + 1.)
        
    if perturbation == True:
        # I replaced the p[25] with inputs_list.
        # inputs_list = Fe_replen_ss = 0 in the ss and
        # inputs_list = MULTIPLE_INPUTS during the experiment
        O2 = O2_INPUT*sigm_time + O2_hypoxia*sigm2_time
    else:
        O2 = O2_INPUT
               
    # Differential equations
    # ddt{Fe}
    ydot[0] = p[4]*Fe_GROWTH*TfR1 - p[2]*Fe - speed_kasoFeFt +\
        (p[14] + p[5])*FtFe + speed_kFtdegFesto
    # ddt{Ft}
    if IRP1 < 1e-30:
        sigp_IRP1_Ft = 0.
        sigp_IRP1_TfR1 = 0.
        sigp_IRP1_HIF = 0.
    else:
        sigp_IRP1_Ft = 1. / ((p[16]/IRP1)**n + 1.)
        sigp_IRP1_TfR1 = 1. / ((p[17]/IRP1)**n + 1.)
        sigp_IRP1_HIF = 1. / ((p[29] / IRP1)**n + 1.)

    if HIF < 1e-30:
        sigp_HIF_TfR1 = 0.
    else:
        sigp_HIF_TfR1 = 1. / ((p[27] / HIF)**n + 1.)

    # ddt{Ft}
    ydot[1] = p[6] - p[31]*sigp_IRP1_Ft - p[5]*Ft - speed_kasoFeFt +\
        speed_kdispkcat_FtFe
        
    # ddt{FtFe}
    ydot[2] = speed_kasoFeFt - speed_kdispkcat_FtFe - p[5]*FtFe
    
    # ddt{Fe_sto}
    ydot[3] = p[13]*FtFe - speed_kFtdegFesto
    
    # ddt{IRP1}
    sigm_ACO1_IRP1 = 1. / ((Fe / p[15])**n + 1.)
    if Fe < 1e-30:
        sigp_IRP1_ACO1 = 0.
    else:
        sigp_IRP1_ACO1 = 1. / ((p[32] / Fe)**n + 1.)
    ydot[4] = p[33] + p[1]*sigm_ACO1_IRP1*(p[0] - IRP1) -\
        p[7]*sigp_IRP1_ACO1*IRP1 - p[8]*IRP1

    # ddt{TfR1}
    ydot[5] = p[11] + p[9]*sigp_IRP1_TfR1 + p[28]* sigp_HIF_TfR1 - p[10]*TfR1

    # ddt{HIF}
    HIF_PHD = 0.5* (HIF-PHD_tot-p[25]+sqrt((p[25]+PHD_tot-HIF)**2 + 4*p[25]*HIF))
    nu_hif_hy = p[24]*PHD_tot*(1/(p[25] + HIF_PHD + PHD_tot))*(O2/(p[26]+O2))
   
    ydot[6] = p[22] - HIF*(p[23] + nu_hif_hy)-p[30]*sigp_IRP1_HIF
    if t >=10000 and t <= 30000.0:    
        JacobianM = [[-(p[2]+p[12]*Ft), 0, p[14]+p[5], p[5], 0, p[4]* Fe_GROWTH, 0],
                     [0, -(p[5]+p[12]*Fe), p[14]+p[13],0, -p[31]*sigp_IRP1_Ft,0,0 ],
                     [0, p[12]*Ft,-p[5]-p[14]-p[13],0,0,0,0],
                     [0,0, p[13], -p[5],0,0,0],
                     [p[1]*sigm_ACO1_IRP1*p[0],0,0,0,-p[1]*sigm_ACO1_IRP1-p[7]*sigp_IRP1_ACO1*IRP1-p[8],0,0],
                     [0,0,0,0,p[9]*sigp_IRP1_TfR1,-p[10],p[28]* sigp_HIF_TfR1],
                     [0, 0, 0,0,p[30]*sigp_IRP1_HIF,0,p[23] + nu_hif_hy]]
        result =la.eig(JacobianM)
        #print "TIME=",t
        result0= result[0].real
        list0 = [result0.tolist()]
        lambda_positive = True
        for jj in list0:
            if jj >= 0.0:
                lambda_positive = True
                break
        if lambda_positive == False:       
            list1 =result[1].tolist()
            list_eigv= list0+list1
            list_eigvT = zip(*list_eigv)
            name_col = ['eigenvectors_' + str(jj) for jj in range(1, len(VAR_NAMES) + 1)]
            name_col = ["eigenvalues" ]+ name_col
            list_eigen = pd.DataFrame(list_eigvT, columns=name_col)
            
            #print list_eigen
            list_R_stifness = []
            name_index = []
            for i in range(0, len(result0)):
                for j in range(i+1, len(result0)):
                    Ratio_stifness = abs(result0[i])/abs(result0[j])
                    list_R_stifness.append(Ratio_stifness)
                    name_index.append("R_stifness_"+str(i)+"_"+str(j))
                    #print "R_stifness_"+str(i)+"_"+str(j)+"=", Ratio_stifness
            #print "List R_stifness", list_R_stifness
            list_R_stiff=pd.DataFrame(list_R_stifness, index= name_index, columns=["Value"])
            #print list_R_stiff
        
            for kk in list_R_stifness:
                if kk > 0:
                    fileStiff.write('It is stiff at time ' + str(t) +'\n')
                    fileStiff.write(pd.DataFrame(list_eigvT, columns=name_col).to_string()+'\n')
                    fileStiff.write(pd.DataFrame(list_R_stifness, index= name_index, columns=["Value"]).to_string() +'\n')

                    
    #********************************
    if t==0:
        min_ydot = ydot 
        min_ydot = ydot 
    else:
        for i in range(0,7):
            if abs(ydot[i]) < abs(min_ydot[i]):
               min_ydot[i] = abs(ydot[i])
            if abs(ydot[i]) > abs(max_ydot[i]):
               max_ydot[i] = abs(ydot[i])
   
    return ydot
# -----------------------------------------------------------------------------


def derivative(f,a,method='central',h=INC):
    '''Compute the difference formula for f'(a) with step size h.

    Parameters
    ----------
    f : function
        Vectorized function of one variable
    a : number
        Compute derivative at x = a
    method : string
        Difference formula: 'forward', 'backward' or 'central'
    h : number
        Step size in difference formula

    Returns
    -------
    float
        Difference formula:
            central: f(a+h) - f(a-h))/2h
            forward: f(a+h) - f(a))/h
            backward: f(a) - f(a-h))/h            
    '''
    if method == 'central':
        return (f(a + h) - f(a - h))/(2*h)
    elif method == 'forward':
        return (f(a + h) - f(a))/h
    elif method == 'backward':
        return (f(a) - f(a - h))/h
    else:
        raise ValueError("Method must be 'central', 'forward' or 'backward'.")

# -----------------------------------------------------------------------------


"""Look for for the first index of a point in a table.
    TODO: Note this function is specific to time table in hours format
    INPUTS:
            point   :: tpoint
            table   :: t_range
    OUTPUT
            index   :: index
"""


def compute_index(tpoint, t_range):

    item = np.where(t_range == tpoint*3600)
    # -> the index value is item[0][0];
    # the corresponding time is thus t_range[item[0][0]]
    # -> the concentration C_i at this time is t_course[index][i],
    # with index = item[0][0].

    try:
        index = item[0][0]
    except IndexError:
        print(('\nERROR in compute_index():  the time point ', tpoint,
              'h does not exist in t_range\n'))
        raise             # To output the stack trace (and exit)
    return index
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------


"""List of the constants defined by the differential equations.
These constant should always be defined after the definition of the
differential equations system. The constants stated as [MANDATORY] are used in
the method module (simulate.py), and should be defined.

START              :: First time instant (in general 0.0)
INC                :: Time step increment
END_SS             :: End of the time interval for the steady state search
END_MAIN           :: End of the time interval for the main simulation
SS_SIMU_ARGS       :: [MANDATORY] informations passed to simulate.py and
                                  containing the steady state odes and time
                                  range
                            --> Pair (Differental equation funtion, time range)
MAIN_SIMU_ARGS     :: [MANDATORY] informations passed to simulate.py and
                                  containing the perturbed odes and time range
                            --> Pair (Differental equation funtion, time range)
INDEX_END          :: index in the time range of the simulation end (53h taking
                      into account the additional 5h of cutoff)
INDEX_DEPRIV_END   :: index in the time range of the deprivation end (25h
                      taking into account the additional 5h of cutoff)
"""

# TODO initialize

SS_SIMU_ARGS = (diff_eqs, np.arange(START, END_SS + INC, INC))
MAIN_SIMU_ARGS = (diff_eqs, np.arange(START, END_MAIN + INC, INC))
# IP 16/04/20 INDEX_START 
INDEX_START = compute_index(0.0, MAIN_SIMU_ARGS[1])
INDEX_END = compute_index(53., MAIN_SIMU_ARGS[1])
INDEX_DEPRIV_END = compute_index(25., MAIN_SIMU_ARGS[1])
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
"""Add or replace parameter values at the end of the parameters set.
    The parameters in the list PARAM_ALGEB_DEFS are defined as functions of the
    other parameters

    In other words the writer of this function should ensure that each
    parameter appearing in PARAM_ALGEB_DEFS is assigned a value, i.e. if we
    have 'k37' in PARAM_ALGEB_DEFS we should have an assignment
    k37 = <expression> in the function code.

    INPUTS:
            Parameters set           :: params
            In-place or copy mode   :: copy

    OUTPUTS:
            Completed/corrected parameters set   :: params
            Boolean checking errors             :: err
            Print associated to errors          :: err_print
"""


def add_params(params, copy=False):

    if copy == True:
        params = params.copy()

    err = 0
    err_print = 'No error'

    # Set k_IRP_Ft inferior but close to k_Ft_prod
    try:
        if DEBUG > 2:
            print('Adding k_IRP_Ft ')
            print('before: k_IRP_Ft = ',  params.k_IRP_Ft)
        # so that  k_Ft_prod - k_IRP_Ft  is positive and small \
        # => production of Ft goes down
        params.k_IRP_Ft = 0.975 * params.k_Ft_prod

        if DEBUG > 2:
            print('after: k_IRP_Ft = ',  params.k_IRP_Ft)
    except:
        err = 1
        err_print = 'Error in adding additional parameter k_IRP_Ft'
        params = []
        return params, err, err_print

    # New constraint ensuring equality of thresholds for IRP <-> Aco1 reactions
    try:
        if DEBUG > 2:
            print('Adding theta_IRP1_Aco1 ')
            print('before: theta_IRP1_Aco1 = ', params.theta_IRP1_Aco1)
        params.theta_IRP1_Aco1 = params.theta_Aco1_IRP1
        if DEBUG > 2:
            print('after: theta_IRP1_Aco1 = ', params.theta_IRP1_Aco1)
    except:
        err = 1
        err_print = 'Error in adding additional parameter theta_IRP1_Aco1'
        params = []
        return params, err, err_print

    # New equation determining k_IRP_prod using d(ACO1+IRP1)/dt = 0
    # and k_ACO1_deg = k_IRP1_deg
    try:
        if DEBUG > 2:
            print('Adding k_IRP_prod ')
            print('before: k_IRP_prod = ', params.k_IRP_prod)
        params.k_IRP_prod = params.C*params.k_IRP1_deg
        if DEBUG > 2:
            print('after: k_IRP_prod = ', params.k_IRP_prod)
    except:
        err = 1
        err_print = 'Error in adding additional parameter k_IRP_prod'
        params = []
        return params, err, err_print

    if DEBUG > 2:
        print(('After adding algebraically defined parameteres,\n params = ',
              params))

    # if len(params) == NB_DIM_TOT - 2:
    #     params += [k_IRP_Ft, theta_IRP1_Aco1]
    # elif len(params) == NB_DIM_TOT:
    #    params[NB_DIM_TOT-2:NB_DIM_TOT] = [k_IRP_Ft, theta_IRP1_Aco1]
    # else:
    #     err = 1
    #     err_print = 'Error in adding additional parameter.\n>>>> len(params)\
        # = ', len(params)
    #     err_print = 'NB_DIM_TOT = ', NB_DIM_TOT

    return params, err, err_print
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------


"""Check static constraint on the parameters. Can modify the parameters set if
    needed.

    INPUTS:
            Parameters set           :: params
            In-place or copy mode   :: copy

    OUTPUTS:
        Completed/corrected parameters set                         :: params
        Boolean checking constraints satisfiability               :: sat
        Boolean checking errors                                   :: err
        Print associated to errors and unsatisfiability           :: err_print
"""

PARAMS_CONSTRAINT_BOX[18][1]
def check_static_constraints_params(params, copy=True):
    if copy == True:
        params = params.copy()
#    print 'lunghezza params', len(params)
#    params.remove(params[10])
#    print len(params)
    err = 0
    err_print = 'No error'
    sat = True              # No hard constraint
    try:
        for i in range(NB_SAMPLED_DIM):
            if DEBUG > 2:
                print('param [', i, '] ; ', PARAM_NAMES[i], ' = ', params[i])
                print('Bounds are :  ', PARAMS_CONSTRAINT_BOX[i])
            if params[i] < PARAMS_CONSTRAINT_BOX[i][0] or params[i]\
                    > PARAMS_CONSTRAINT_BOX[i][1]:
                sat = False
                err_print = '{}  = {} :: ouside of bounds: [{}, {}]'
                err_print = err_print.format(
                                             PARAM_NAMES[i],
                                             params[i],
                                             PARAMS_CONSTRAINT_BOX[i][0],
                                             PARAMS_CONSTRAINT_BOX[i][1]
                                             )
                return params, sat, err, err_print
    except:
        err = 1
        sat = False
        err_print = 'Error in checking parameters [box constraints check].\n'
        params = []
        return params, sat, err, err_print
   
    return params, sat, err, err_print
# -----------------------------------------------------------------------------


def algebraic_steady_state(params):
    err = 0
    err_print = 'No error'
    # The code of this function should be generated automatically
    # by SAGE in a separate file
    # Note: it is e value at t=0, when iron is abundant
# This var is used, not assigned to.
# So the global statement is unnecessary here.
    # The parameters are used, but not assigned to  =>  the global statement
    # is not needed for them
    # To be completed if this function is put in the main python file :
    # C = check_with_param_sets.C;
    # k_Aco1_IRP1 = check_with_param_sets.k_Aco1_IRP1
    # k_FPN1a_deg = check_with_param_sets.k_FPN1a_deg;
    # k_FPN1a_prod = check_with_one_exp.k_FPN1a_prod
    # ATTENTION : k_Fe_IRP est obsolete  (nous avons maintenant k_Aco1_IRP1 et
    # k_IRP1_Aco1)
    # -----------------------------------------------------------------------
    # Equations provided by SAGE  (voir define_CI.m/set_to_steady.m de
    # simplif4_5B)
    try:
        C = params[0]
        k_Aco1_IRP1 = params[1]
        k_Fe_cons = params[2]
        k_Fe_export = params[3]
        k_Fe_input = params[4]
        k_Ft_deg = params[5]
        k_Ft_prod = params[6]
        k_IRP1_Aco1 = params[7]
        k_IRP1_deg = params[8]
        k_IRP_TfR1 = params[9]
        k_TfR1_deg = params[10]
        k_TfR1_prod = params[11]
        k_aso = params[12]
        k_cat = params[13]
        k_dis = params[14]
        theta_Aco1_IRP1 = params[15]
        theta_IRP_Ft = params[16]
        theta_IRP_TfR1 = params[17]
        k_FPN1a_deg = params[18]
        k_FPN1a_prod = params[19]
        k_IRP_FPN1a = params[20]
        theta_IRP_FPN1a = params[21]
        k_hif_prod         = params[22]
        k_hif_deg          = params[23]
        k_pdh_cat          = params[24]
        k_phd_dis          = params[25]
        k_phd_M            = params[26]
        theta_HIF_TfR1     = params[27]
        k_HIF_TfR1         = params[28]
        theta_IRP_HIF      = params[29]
        k_IRP_HIF          = params[30]
        k_IRP_Ft = params[31]
        theta_IRP1_Aco1 = params[32]
        k_IRP_prod = params[33]
    except:
        err = 1
        err_print = 'model. algebraic_steady_state():\n         \
                        >>> Error in handling parameter list.\n >>> \
                        Length params list =  {}'
        err_print = err_print.format(len(params))
        if DEBUG > 2:  # Anormal behaviour
            print(err_print)
            assert(0)
        return [], err, err_print
    try:
        # Compacted version
        # IP 05/03/19 I introduce A0, A1,A3, fat1, fat2,fat2a, fat3 fo better
        # read the formules.
        A0 = k_FPN1a_deg*k_Fe_cons
        A1 = k_Ft_prod*k_TfR1_deg
        fat1 = Fe_GROWTH*k_FPN1a_deg*k_Fe_input*k_Ft_prod*k_TfR1_prod*k_aso
        fat2 = Fe_GROWTH*k_FPN1a_deg*k_Fe_input*k_Ft_deg**2*k_TfR1_prod*k_aso
        fat2a = Fe_GROWTH*k_FPN1a_deg*k_Fe_input*k_Ft_deg*k_TfR1_prod*k_aso
        fat3 = k_FPN1a_deg*k_Fe_cons
        A3 = k_Ft_deg + k_cat + k_dis
        Fe_0 = Fe_GROWTH*k_FPN1a_deg*k_Fe_input*k_TfR1_prod/(fat3*k_TfR1_deg)

        Ft_0 = (A0*A1*A3)/(fat2a
                           + fat3*k_Ft_deg**2*k_TfR1_deg
                           + fat3*k_Ft_deg*k_TfR1_deg*k_cat
                           + fat3*k_Ft_deg*k_TfR1_deg*k_dis)

        FtFe_0 = fat1/(fat2a
                       + fat3*k_Ft_deg**2*k_TfR1_deg
                       + fat3*k_Ft_deg*k_TfR1_deg*k_cat
                       + fat3*k_Ft_deg*k_TfR1_deg*k_dis)

        Fe_sto_0 = fat1*k_cat/(fat2
                               + fat3*k_Ft_deg**3*k_TfR1_deg
                               + fat3*k_Ft_deg**2*k_TfR1_deg*k_cat
                               + fat3*k_Ft_deg**2*k_TfR1_deg*k_dis)

        IRP1_0 = k_IRP_prod/(k_IRP1_Aco1 + k_IRP1_deg)

        # FPN1a_0 = k_FPN1a_prod/k_FPN1a_deg
        nu_hif_hy = k_pdh_cat*PHD_tot*(1/(k_phd_dis + HIF_PHD + PHD_tot))*(O2_INPUT/(k_phd_M+O2_INPUT))
        HIF_0 = k_hif_prod / (k_hif_deg + nu_hif_hy)

        TfR1_0 = k_TfR1_prod/k_TfR1_deg
    except:
        err = 1
        err_print = 'model.algebraic_steady_state():\n                    >>> \
            Error in Computing initial state.'
        if DEBUG > 2:                           # Anormal behaviour
            print(err_print)
            assert(0)
        return [], err, err_print

    if DEBUG > 1:
        print(('\ninitial values of [Fe, Ft, FtFe, Fe_sto, IRP1, FPN1a, TfR1]'
              ' =', [Fe_0, Ft_0, FtFe_0, Fe_sto_0, IRP1_0, TfR1_0, HIF_0], '\n'))
    init_vars = [Fe_0, Ft_0, FtFe_0, Fe_sto_0, IRP1_0, TfR1_0, HIF_0]

    return init_vars, err, err_print
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------


"""Check constraints on the steady state initial condition.

    INPUTS:
            Initial conditions  ::  init_vars
            Parameters set       ::  params

    OUTPUT:
            Boolean checking constraints satisfiability     :: sat
            Boolean checking errors                         :: err
            Print associated to errors and unsatifiability  :: err_print
"""


def check_initial_state_constraint(init_vars, params):

    err = 0
    err_print = 'No error\n'
    sat = True

    ivar = 0
    for var in init_vars:
        if var < VARS_BOUNDS[ivar][0] or var > VARS_BOUNDS[ivar][1]:
            sat = False
            err_print += 'The variable \'{}\' = {} is outside of bounds {}'\
                .format(VAR_NAMES[ivar], var, VARS_BOUNDS[ivar])
            return sat, err, err_print
        ivar += 1

    Fe, Ft, FtFe, Fe_sto, IRP1, TfR1, HIF = init_vars[0], init_vars[1],\
        init_vars[2], init_vars[3], init_vars[4], init_vars[5], init_vars[6]

    if not Fe > params.theta_Aco1_IRP1:
        err_print += 'Unsat: Fe > theta_Aco1_IRP1 ---> skipping...\n'
        sat = False

    if not Fe > params.theta_IRP1_Aco1:
        err_print += 'Unsat: Fe > theta_IRP1_Aco1 ---> skipping...\n'
        sat = False

    # Constraints on Initial Ft concentrations (Ft[0])
    if Ft < FT_INIT_CONSTRAINT[0]:  # Relaxed constraint of 5nM instead of 50nM
        err_print += 'Unsat: Ft(0) < 5e-09 mol/L ---> skipping...\n'
        sat = False

    if Ft > FT_INIT_CONSTRAINT[1]:
        err_print += 'Unsat: Ft(0) > 1.0e-06 mol/L ---> skipping...\n'
        sat = False
    # Constraints on Initial HIF concentrations (HIF[0])
    if HIF < HIF_INIT_CONSTRAINT[0]:
        err_print += 'Unsat: HIF(0) <'+ str(HIF_INIT_CONSTRAINT[0])+ ' mol/L ---> skipping...\n'
        sat = False

    if HIF > HIF_INIT_CONSTRAINT[1]:
        err_print += 'Unsat: HIF(0) >'+ str(HIF_INIT_CONSTRAINT[1])+ ' mol/L ---> skipping...\n'
        sat = False
    # For IRP1:
    #     sigp(IRP1,theta_IRP_Ft,n) should be close to 0, and
    #     sigp(IRP1,theta_IRP_FPN1a,n) should be close to 0
    if not IRP1 < params.theta_IRP_Ft:
        err_print += 'Unsat: IRP1 < theta_IRP_Ft ---> skipping...\n'
        sat = False

    # Similar constraint for TfR1
    if not IRP1 < params.theta_IRP_TfR1:
        err_print += 'Unsat: IRP1 < theta_IRP_TfR1 ---> skipping...\n'
        sat = False

    # *** Check the other constraints bearing on the steady state
    #     ( Then, ideally:  check_steady_constraints(PRLog); )
    #     ( See  .../breach/check_steady_constraints.m )

# TODO This constraint was never check (relaxed or forsaken)
        # sat = False

    if not params.C > IRP1:
        err_print += 'Unsat: C > IRP1'
        sat = False
    # TODO also check Aco1 = C-IRP1 > IRP1 at the initial condition
    stored_iron = Fe_sto/(Ft+FtFe)
    # 50 <= Fe_sto/(Ft+FtFe) <= 2500
    if not stored_iron <= Fe_By_Ft_CONSTRAINT[1] and not \
            stored_iron >= Fe_By_Ft_CONSTRAINT[0]:
        err_print += 'Unsat: 50 <= Fe_sto/(Ft+FtFe) <= 2500\n'
        sat = False
    return sat, err, err_print
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------


"""For ONE input (here one value of Fe_replen) check constraints on the
    results of the associated simulation
    INPUTS:
            Simulation results  :: simulation
                -> Tuple: (traces_variables, time_points, initial_conditions)
            Current input value ::  current_input

    OUTPUTS:
        Boolean checking the satisfiability of the constraints  :: sat
        Boolean checking if errors occured                      :: err
        Print associated to errors and unsatifiability          :: err_print
"""


def check_absolute_constraints(simulation, current_input):  # TODO

    sat = True
    err = 0
    err_print = ''

    # Reminder: simulation = (traces_variables, time_points,\
    # initial_conditions)
    traces = simulation[0]
    time = simulation[1]
#    init_cond = simulation[2]  # IP never used

# INDEX_END is defined at the module initialization.
    Fe_end = traces[INDEX_END][0]
    # Currently INDEX_END correspond to t = 53h as t_cutoff ~ 5h
    # Ft[end of re-entry exp]
    Ft_end = traces[INDEX_END][1]
    # Ft:Fe[end of re-entry exp]
    FtFe_end = traces[INDEX_END][2]
    # Fe-sto[end of re-entry exp]
    Festo_end = traces[INDEX_END][3]
    # (Ft+Ft:Fe)[end of re-entry exp]
    # IRP[end of re-entry exp]
    # (Fe + Ft:Fe + Fe-sto)[end of re-entry exp]
    Fetot_end = Fe_end + FtFe_end + Festo_end
    # FPN1a[end of re-entry exp]
    TfR1_end = traces[INDEX_END][5]
    # IP 04/10/19 - Need to add absolute constraint for HIF
    HIF_start = traces[INDEX_START][6]
    HIF_end = traces[INDEX_END][6]
    # TfR1[end of re-entry exp]
    # INDEX_DEPRIV_END is defined at the module initialization
    Fefree_depriv_end = traces[INDEX_DEPRIV_END][0]
    # Currently INDEX_DEPRIV_END correspond to t = 25h

    if DEBUG > 2:
        print(('## End time (', INDEX_END, ') solutions = \n',
              traces[INDEX_END], '\n'))
        print('Replen = ', current_input)
        print('Fe_end =', traces[INDEX_END][0])
        print('Ft_end   = ', traces[INDEX_END][1])
        print('FtFe_end = ', traces[INDEX_END][2])
        print('Festo_end = ', traces[INDEX_END][3])
        print('Fttot_end= ', Ft_end + FtFe_end)
        print('IRP_end =  ', traces[INDEX_END][4])
        print('Fetot_end = ', Fe_end + FtFe_end + Festo_end)
        print('TfR1_end  =  ', traces[INDEX_END][5])
        print('HIF_end = ', traces[INDEX_END][6])
        print('Fefree_depriv_end = ', traces[INDEX_DEPRIV_END][0], '\n')


    # ********* CONSTRAINTS CHECKS *************
    # ********* Check behavioral constraints ***
    # alw (Fe[t] < 2.0e-6)   ref 39 du draft, Epsztejn, , Cabantchik, 1997
    # Max bound on Fe for biological viability
    list_invalid_Fe = [
                       traces[t_index][0]
                       for t_index, t_value in enumerate(time)
                       if traces[t_index][0] > MAX_Fe_CONCENTRATION
                       ]
    if len(list_invalid_Fe) > 0:
        err_print += 'Unsat:  alw (Fe[t] < 2.0e-5)\n'
        err_print += '    Values found: {}\n'.format(list_invalid_Fe)
        sat = False
        
    list_invalid_HIF = [
            traces[t_index][0]
            for t_index, t_value in enumerate(time)
            if traces[t_index][6] > MAX_HIF_CONCENTRATION
            ]
    if len(list_invalid_HIF) > 0:
        err_print += 'Unsat:  alw (HIF[t] < 2.0e-6)\n'
        err_print += '    Values found: {}\n'.format(list_invalid_HIF)
        sat = False
#    list_invalid_FtFe = [
#            traces[t_index][0]
#            for t_index, t_value in enumerate(time)
#            if traces[t_index][6] > MAX_FtFe_CONCENTRATION
#            ]
#    if len(list_invalid_FtFe) > 0:
#        err_print += 'Unsat:  alw (HIF[t] < 2.0e-6)\n'
#        err_print += '    Values found: {}\n'.format(list_invalid_FtFe)
#        sat = False
    if abs(HIF_end/HIF_start)< 2.5:
        err_print += 'Unsat:  HIF[t] costant line \n'
        sat = False
        
#    if current_input >= 2.5e-07 and (
#                                     Fetot_end < CSTRT_Fetot_HIGH_END[0]
#                                     or Fetot_end > CSTRT_Fetot_HIGH_END[1]
#                                     ):
#        err_print += 'Unsat Fetot[48h]_Tf125 = {} :  > 2000e-09 mol/L ' +\
#            'or < 5e-09 mol/L  ---> skipping\n'.format(Fetot_end)
#        sat = False
#
#    if current_input == 0.0 and (Fetot_end > CSTRT_Fetot_LOW_END):
#        err_print += 'Unsat Fetot[48h]_Tf0 =  {} : > 15e-09 mol/L ' + \
#            '---> skipping\n'.format(Fetot_end)
#        sat = False

    if Fefree_depriv_end < MIN_FREE_IRON:
        err_print += 'Unsat Fefree[24h] =  {} : < 1e-11 mol/L ' + \
            '---> skipping\n'.format(Fefree_depriv_end)
        sat = False

    if sat == True and err == 0:
        err_print = 'No error'
    return sat, err, err_print
# -----------------------------------------------------------------------------


# -----------------------------------------------------------------------------
# TODO: ERROR are not well handled


"""Compute the relative value from the results of simulations with various
    inputs Here we compute the ratio Var_replen[END] / Var_no-replen[END], with
    Var bing Ft and TfR1
    TODO: Note there is 3 ratios instead of 2 to be coherent with old filter.py
    code
    INPUT:
            List of simulation results for each inputs  :: simulations_list

    OUTPUTS:
            Computed relative value (ratio) :: relative_sols
                    -> List of tuples [ (ratio1_input1, ratio2_input1, ...),
                                        (ratio1_input2, ratio2_input2, ...),
                                        ... ]
            Boolean checking if errors occured      :: err
            Print associated to errors              :: err_print
"""


def compute_relative_values(simulations_list):
    err = 0
    err_print = 'No error' # TODO Add div by zero error handling ?
    relative_sols = []
    i = 0

    # No rentry simulation
    simulation_LOW = simulations_list[-1]
    traces_LOW = simulation_LOW[0]
#    time_LOW = simulation_LOW[1]  # IP never used
#    init_cond_LOW = simulation_LOW[2]  # IP never used
    # Ft[end of re-entry exp]
    Ft_end_LOW = traces_LOW[INDEX_END][1]
    # Ft:Fe[end of re-entry exp]
    FtFe_end_LOW = traces_LOW[INDEX_END][2]
    # (Ft+Ft:Fe)[end of re-entry exp]
    Fttot_end_LOW = Ft_end_LOW + FtFe_end_LOW
    # TfR1[end of re-entry exp]
    TfR1_end_LOW = traces_LOW[INDEX_END][5]

    for simulation in simulations_list[:len(MULTIPLE_INPUTS)-1]:
        # Reminder: simulation = (traces_variables, time_points,
        # initial_conditions)
        traces = simulation[0]
#        time = simulation[1]  # IP never used
#         init_cond = simulation[2]  # IP never used
        # Ft[end of re-entry exp]
        Ft_end = traces[INDEX_END][1]
        # Ft:Fe[end of re-entry exp]
        FtFe_end = traces[INDEX_END][2]
        # (Ft+Ft:Fe)[end of re-entry exp]
        Fttot_end = Ft_end + FtFe_end
        # TfR1[end of re-entry exp]
        TfR1_end = traces[INDEX_END][5]

        # TfR1 Ratio
        if DEBUG > 1:
            print(('TfR1_end[{}]/TfR1_end[{}] = {}/{} \n'
                  .format(i, len(MULTIPLE_INPUTS)-1, TfR1_end, TfR1_end_LOW)))
        try:
            ratio_TfR1 = TfR1_end/TfR1_end_LOW
        except:
            err = 1
            err_print = 'Numerical error in ratio computation'
            return [], err, err_print
        if isnan(ratio_TfR1) or isinf(ratio_TfR1):
            err = 1
            err_print = 'Numerical error in ratio computation'
            return [], err, err_print

        # Ft Ratio
        if DEBUG > 1:
            print(('Fttot_end[{}] = {} / Fttot_end[{}] = {} '
                  .format(i, Fttot_end, len(MULTIPLE_INPUTS)-1, Fttot_end_LOW)))
        try:
            ratio_Ft = Fttot_end/Fttot_end_LOW
        except:
            err = 1
            err_print = 'Numerical error in ratio computation'
            return [], err, err_print
        if isnan(ratio_Ft) or isinf(ratio_Ft):
            err = 1
            err_print = 'Numerical error in ratio computation'
            return [], err, err_print
        i += 1
        # Storing ratios
        # TODO: Note the -1 for the used Fe-sto ratios...
        # Temporary....
        relative_sols.append((-1.0, ratio_TfR1, ratio_Ft))

    return relative_sols, err, err_print
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------


"""Check constraint on relative results.
    Prototype function: EMPTY
"""


def check_relative_constraints(relative_result_list):
    sat, err, err_print = True, 0, 'No error'
    if DEBUG > 0:
        print('NO RELATIVE CONSTRAINTS CHECK')
    return sat, err, err_print
# -----------------------------------------------------------------------------


# -----------------------------------------------------------------------------
# TODO: ERROR HANDLING NOT DONE


"""Create a string output which summarize the results.
    The format is the one read by filter.py

    INPUTS:
            Simulations solutions for one parameters set :: sol
    OUTPUT:
        String containing the information to write in file  :: formated_string
        Boolean checking if an error occured                :: err
        Print associated to errors                          :: err_print
"""


def format_solution(sol):

    # Reminder
    # sol = (params, relative_sols, all_simus, init_vars)  :: Tuple
    # relative_sols = [(-1.0, ratio_TfR1, ratio_Ft), ...]  :: list of Tuples
    # (NB INPUTS)
    # all_simus = [ (trace, t_range, init_vars),...]       :: List of Tuples
    # trace = Arrays.shape(lines=nb_steps, cols=nb_vars)   :: Array
    # TODO error handling
    err = 0
    err_print = 'No error handling written'

    # FORMAT FOR PRINTNG IN FILE #####
#    nb_params = len(sol[0])  # IP never used
    nb_init_vars = len(sol[3])
    fmt_sols_1 = '\n' + '{} '*(NB_DIM_TOT) + '\n'  # Saving parameters
    fmt_sols_1 = fmt_sols_1.format(*sol[0])

    fmt_sols_2 = '    ' + '{} '*(nb_init_vars) + '\n'  # Saving initial values
    fmt_sols_2 = fmt_sols_2.format(*sol[3])

    fmt_sols_3 = '    e{} ' + '{} '*(NB_RATIOS) + '\n'
    fmt_ratio = ''
    for i in range(len(MULTIPLE_INPUTS)-1):
        if i == 0:
            fmt_ratio = fmt_sols_3.format((i+1), *sol[1][i])
        else:
            fmt_ratio = fmt_ratio + fmt_sols_3.format((i+1), *sol[1][i])

# TODO Output cost function in case of optimization, currently not used
    fmt_cost = 'J= ' + '{}\n'.format(-1.)

    fmt_Fetot = 'Fetot(250,25,10,0)= '
    all_simus = sol[2]
    for i in range(len(MULTIPLE_INPUTS)):
        traces = all_simus[i][0]
        # Currently INDEX_END correspond to t = 53h as t_cutoff ~ 5h
        Fe = traces[INDEX_END][0]
        # Ft:Fe[end of re-entry exp]
        FtFe = traces[INDEX_END][2]
        Festo = traces[INDEX_END][3]
        # (Fe + Ft:Fe + Fe-sto)[end of re-entry exp]
        Fetot = Fe + FtFe + Festo
        s_Fe_tot = '{} '.format(Fetot)
        fmt_Fetot = fmt_Fetot + s_Fe_tot
    fmt_Fetot = fmt_Fetot + '\n'

    formated_string = fmt_sols_1 + fmt_sols_2 + fmt_ratio + fmt_cost +\
        fmt_Fetot

    return formated_string, err, err_print
# -----------------------------------------------------------------------------


# -----------------------------------------------------------------------------
# TODO: Error handling not done
"""Console print additional results which are model specific.
    Here it plots the ratio list.

    INPUT:
            The list of all the solutions   :: list_sols

    OUTPUTS:
            Boolean checking if an error occured    :: err
            Print associated to errors              :: err_print
"""


def print_additional_results(list_sols):

    err = 0
    err_print = 'No error handling ...'  # TODO

    list_rel_sols = [sol[1] for sol in list_sols]
    print('\nSolution relative constraints bounds: ')
    for i in range(len(MULTIPLE_INPUTS)-1):
        TfR1_ratios = [rel_sols[i][1] for rel_sols in list_rel_sols]
        Ft_ratios = [rel_sols[i][2] for rel_sols in list_rel_sols]
        bl_rTfR1 = np.amin(np.array(TfR1_ratios))
        bu_rTfR1 = np.amax(np.array(TfR1_ratios))
        bl_rFt = np.amin(np.array(Ft_ratios))
        bu_rFt = np.amax(np.array(Ft_ratios))
        print(('Experiment {} ratio TfR1:  {}    {}'
              .format(i+1, bl_rTfR1, bu_rTfR1)))
        print(('Experiment {} ratio Ft:    {}    {}\n'
              .format(i+1, bl_rFt, bu_rFt)))
    return err, err_print
