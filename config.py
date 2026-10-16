# -*- coding: utf-8 -*-
"""
Created on Mon Mar 30 14:15:00 2020

@author: irina
"""


""" config.py will be stored all the constant variable
"""
# ________________________________Simulate_____________________________________
DEBUG_S = 2
MAX_ITERATION_SS_IMPROVEMENT = 20  # 20  # default 10
SIMU_RTOL = 1e-8    # default 1e-10 # low 10e-8
SIMU_ATOL = 1e-10    # default 1e-13 # low 10e-10
EPSILON = 8e-6  # 1e-8  # 8e-6
# IP 11/03/19
# Nb steady points to consider that the system is in a steady state
# It means that this 63 points define almost a costant line because the
# system does not change its state.
NB_STEADY_POINTS = 63
PLOT_MODE = 'traces'  # 'ribbon' or 'traces'
NUMERICAL_SOLVER = 'scipy-odeint-default'
# NUMERICAL_SOLVER= 'scipy-solve_ivp-LSODA'
# NUMERICAL_SOLVER = 'julia'
FORCED_SAMPLING_SCALE = ''  # '' by default or 'linear'
#FORCED_SAMPLING_SCALE = 'linear'
# moved from model to simulate
ADD_SAMPLE_FILE = 'temp.box'
NB_ADD_SAMPLES = 1000
RESAMPLING_SIZE = 0.4  # or 0.4 or 0.2
skip = 1280
magnitude = 2
increment = 1000
PCA_REDUCTION = 'reduced'  # 'reduced' or 'full'
# Cap how many full ODE trajectories are re-simulated / dumped / plotted
# after sampling (keeps post-processing memory bounded).
MAX_TRAJECTORIES_DUMP = 100
MAX_PLOT_CURVES = 20
# _____________________________model.py_______________________________________
NB_PLOT_POINT = 5000  # unused legacy alias; see MAX_PLOT_CURVES
Name_Model = 'A model of the Cellular HIF - Iron Homeostasis Network'
map_VAR_NAMES = {'Fe'     : 1,
                 'Ft'     : 2,
                 'FtFe'   : 3,
                 'Fe_sto' : 4,
                 'IRP1'   : 5,
                 'TfR1'  : 6,
                 'HIF'   : 7
                 }
map_PARAM_NAMES = {'C'              : 1,
                   'k_Aco1_IRP1'    : 2,
                   'k_Fe_cons'      : 3,
                   'k_Fe_export'    : 4,
                   'k_Fe_input'     : 5,
                   'k_Ft_deg'       : 6,
                   'k_Ft_prod'      : 7,
                   'k_IRP1_Aco1'    : 8,
                   'k_IRP1_deg'     : 9,
                   'k_IRP_TfR1'     : 10,
                   'k_TfR1_deg'     : 11,
                   'k_TfR1_prod'    : 12,
                   'k_aso'          : 13,
                   'k_cat'          : 14,
                   'k_dis'          : 15,
                   'theta_Aco1_IRP1': 16,
                   'theta_IRP_Ft'   : 17,
                   'theta_IRP_TfR1' : 18,
                   'k_FPN1a_deg'    : 19,
                   'k_FPN1a_prod'   : 20,
                   'k_IRP_FPN1a'    : 21,
                   'theta_IRP_FPN1a': 22,
                   'k_hif_prod'     : 23,
                   'k_hif_deg'      : 24,
                   'k_pdh_cat'      : 25,
                   'k_phd_dis'      : 26,
                   'k_phd_M'        : 27,
                   'theta_HIF_TfR1' : 28,
                   'k_HIF_TfR1'     : 29,
                   'theta_IRP_HIF'  : 30,
                   'k_IRP_HIF'      : 31
                   }
Name_Model = 'A model of the Cellular HIF - Iron Homeostasis Network'
VAR_NAMES = ['Fe',
             'Ft',
             'FtFe',
             'Fe-sto',
             'IRP1',
             'TfR1',
             'HIF']
PARAM_NAMES = ['C',                # 1
               'k_Aco1_IRP1',      # 2
               'k_Fe_cons',        # 3
               'k_Fe_export',      # 4
               'k_Fe_input',       # 5
               'k_Ft_deg',         # 6
               'k_Ft_prod',        # 7
               'k_IRP1_Aco1',      # 8
               'k_IRP1_deg',       # 9
               'k_IRP_TfR1',       # 10
               'k_TfR1_deg',       # 11
               'k_TfR1_prod',      # 12
               'k_aso',            # 13
               'k_cat',            # 14
               'k_dis',            # 15
               'theta_Aco1_IRP1',  # 16
               'theta_IRP_Ft',     # 17
               'theta_IRP_TfR1',   # 18
               'k_FPN1a_deg',      # 19
               'k_FPN1a_prod',     # 20
               'k_IRP_FPN1a',      # 21
               'theta_IRP_FPN1a',  # 22
               'k_hif_prod',       # 23
               'k_hif_deg',        # 24
               'k_pdh_cat',        # 25
               'k_phd_dis',        # 26
               'k_phd_M',          # 27
               'theta_HIF_TfR1',   # 28
               'k_HIF_TfR1',       # 29
               'theta_IRP_HIF',    # 30
               'k_IRP_HIF',        # 31

               'k_IRP_Ft',         # 32
               'theta_IRP1_Aco1',  # 33
               'k_IRP_prod'        # 34
               ]
PARAM_ALGEB_DEFS = ['k_IRP_Ft',
                    'theta_IRP1_Aco1',
                    'k_IRP_prod']
NB_DIM_TOT =               34
NB_SAMPLED_DIM =           31
# 19 # without FPN1a parameters to sample
# 23 # with FPN1a parameters to sample
DEFAULT_ADD_PARAMS_NAMES = []  # not used
DEFAULT_ADD_PARAMS = [[]]
# [[2.e-6,5.e-16,5.e-16,2.e-6]]  # FPN1a parameteres are fixed by default
# [[]]     # FPN1a parameters are not fixed by default
NB_VARS = 7
Fe_replen_ss = 0.0
# MULTIPLE_INPUTS = [250e-09]
MULTIPLE_INPUTS = [0.0]
NB_INPUTS = len(MULTIPLE_INPUTS)
Fe_GROWTH = 2.4e-06
TIME_CUTOFF = 5.0
O2_INPUT =                200e-06  # IP we need to search for the correct value
O2_hypoxia =              9.9e-6
PHD_tot =                 1.0e-06  # IP we need to search for the correct value
HIF_PHD =                 0.0      # IP we need to search for the correct value
VARS_BOUNDS = [
        (1.0e-10, 2.0e-06),  # Fe
        (5.0e-09, 1.0e-06),   # Ft
        (1.0e-15, 5.0e-05),    # Ft:Fe
        (1.0e-11, 5.0e-05),   # Fe-sto
        (1.0e-10, 5.0e-08),   # IRP1
        (1.0e-08, 3.0e-05),   # TfR1
        (1.0e-10, 5.0e-04)    # HIF
        ]

# IP 13 03 19 PARAMS_CONSTRAINT_BOX acts on the first 22 parameters.
PARAMS_CONSTRAINT_BOX = [
        (1.0e-10, 1.0e-04),  # 1  C                  concentration
        (1.0e-07, 1.0e-01),  # 2  k_Aco1_IRP1        1/time
        (1.0e-09, 1.0e03),   # 3  k_Fe_cons          1/time
        (1.0e-15, 1.0e16),   # 4  k_Fe_export        1/(time*concentration)
        (1.0e01, 1.0e08),    # 5  k_Fe_input         1/time
        (3.0e-08, 4.0e-03),  # 6  k_Ft_deg           1/time
        (1.0e-14, 2.0e-09),  # 7  k_Ft_prod          concentration/time
        (5.0e-08, 5.0e-01),  # 8  k_IRP1_Aco1        concentration/time
        (1.0e-08, 2.0e-03),  # 9  k_IRP1_deg         1/time
        (1.0e-19, 1.0e-07),  # 10 k_IRP_TfR1         concentration/time
        (5.0e-08, 1.0e-02),  # 11 k_TfR1_deg         1/time
        (1.0e-16, 1.0e-05),  # 12 k_TfR1_prod        concentration/time
        (2.0e02, 2.0e11),    # 13 k_aso              1/time
        (1.0e01, 1.0e06),    # 14 k_cat              1/time
        (1.0e00, 1.0e06),    # 15 k_dis              1/time
        (1.0e-11, 1.0e-06),  # 16 theta_Aco1_IRP1    concentration
        (1.0e-11, 1.0e-06),  # 17 theta_IRP_Ft       concentration
        (1.0e-11, 1.0e-06),  # 18 theta_IRP_TfR1     concentration
        (1.0e-11, 1.0e-03),  # 19 k_FPN1a_deg        1/time
        (1.0e-19, 1.0e-05),  # 20 k_FPN1a_prod       concentration/time
        (1.0e-20, 1.0e-05),  # 21 k_IRP_FPN1a        concentration/time
        (1.0e-12, 1.0e-03),  # 22 theta_IRP_FPN1a    concentration
        (1.0e-16, 1.0e-09),   # 23 k_hif_prod         concentration/time
        (1.0e-08, 1.0e+01),   # 24 k_hif_deg          1/time
        (1.0e+00, 1.0e+05),   # 25 k_pdh_cat          1/time
        (1.0e-02, 1.0e+05),   # 26 k_phd_dis          concentration
        (1.0e-06, 1.0e-02),   # 27 k_phd_M            concentration
        (1.0e-11, 1.0e-06),   # 28 theta_HIF_TfR1     concentration
        (1.0e-13, 1.0e-06),   # 29 k_HIF_TfR1         concentration/time
        (1.0e-11, 1.0e-06),   # 30 theta_IRP_HIF      concentration
        (1.0e-13, 1.0e-06)    # 31 k_IRP_HIF          concentration/time
        ]

# Relaxed constraint of 5nM instead of 50nM
FT_INIT_CONSTRAINT = [5e-09, 1.0e-06]  # Duplicate of above
HIF_INIT_CONSTRAINT =  [1e-12, 1e-9]
Fe_By_Ft_CONSTRAINT = [50, 2500]
# IP 09/10/19 I changed the max iron concentration of one order of magnitude
MAX_Fe_CONCENTRATION = 2.0e-05
MAX_HIF_CONCENTRATION = 2.0e-06
MAX_FtFe_CONCENTRATION = 1.5e-12

# IP 09/10/19 I changed the max iron concentration upper limit
# of one order of magnitude
CSTRT_Fetot_HIGH_END = [5e-09, 2.5e-05]
CSTRT_Fetot_LOW_END = 15e-09
MIN_FREE_IRON = 1e-11
NB_RATIOS = 3
DEBUG_M = 1  # TODO set in initializig function
# time settings
START = 0.0
INC = 600.  # step 10 minutes
INC_m = INC/60
END_SS = 50.*3600  # 50 hours
END_SS_m = 50
END_MAIN = 75.*3600  # 53 hours
END_MAIN_m = 75

# _____________________________Filter.py_______________________________________
# ----- Values of the bounds for the constraints on Fe_sto, TfR1 and Ft
# OLD Constraints Fe_sto TODO delete this in code....
constraints_on_Fe_sto = [(-1e30, 1e30),
                         (-1e30, 1e30),
                         (-1e30, 1e30)]

# Validated constraints 07/08/2018
# constraints_on_TfR1 = [(0.4, 0.6),
#                       (0.7, 1.0),
#                       (0.8, 1.0)]
# constraints_on_Ft = [(4.0, 7.0),  #OLD Constraint on Ft (relaxed by JMM)
#                     (2.0, 4.0),
#                     (1.0, 1.5)]
# constraints_on_Ft = [(4.0, 12.0),
#                     (2.0, 4.0),
#                     (1.0, 3.0)]

# Relaxed Constraint
constraints_on_TfR1 = [(0.4, 0.6),
                       (0.7, 1.0),
                       (0.8, 1.0)]
constraints_on_Ft = [(4.0, 12.0),
                     (2.0, 4.0),
                     (1.0, 3.0)]

# Constraint (03/09/2018)
# constraints_on_TfR1 = [(0.4, 0.6),
#                       (0.7, 1.0),
#                       (0.8, 1.0)]
# constraints_on_Ft = [(4.0, 12.0),
#                     (2.0, 4.0),  # Relax
#                     (1.0, 3.0)]
Constraints_on_Ft0 = [50.0e-09, 1.0e-06]
Constraints_on_Fetot = [15e-09, 50e-09, 2.5e-06]
# Constraints_on_Ft0 = [5.0e-08,0.5e-06]
# Constraints_on_Fetot = [5e-09,50e-09,0.5e-06]
upperbound_on_C = [1.0e-05]
