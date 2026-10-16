#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 28 12:05:01 2019

@author: piazzai
"""
from __future__ import print_function

import pandas as pd
import sys
import inspect

import config as config
from config import DEBUG_S as debug

ee = 1e-10


def main(filename=None):
    if filename is None:
        filename = sys.argv[1]
    params_listT = read_params_results(filename)
    print(params_listT, '\n')


def read_file(filename, debug, map_name):
    """
    Read variable/parameter bound files and return sorted (lb, ub, name) tuples.
    """
    try:
        data = pd.read_csv(filename,
                           usecols=['name',
                                    'lower_bounds',
                                    'upper_bounds'],
                           sep=r'\s+',
                           engine='python')
    except IOError:
        func_name = inspect.currentframe().f_code.co_name
        sys.stderr.write(
            '\n ERROR in %s: cannot open file %s\nabort...\n'
            % (func_name, filename))
        sys.exit(1)
    data.columns = data.columns.str.strip()
    N_par = len(map_name)
    newcol = []
    absent_col = list(range(1, N_par + 1))
    name_by_index = dict((v, k) for k, v in map_name.items())
    for index, row in data.iterrows():
        newcol = newcol + [map_name[row['name']]]
        absent_col.remove(map_name[row['name']])
    for ii in absent_col:
        data = pd.concat(
            [data, pd.DataFrame([{
                'name': name_by_index[ii],
                'lower_bounds': ee,
                'upper_bounds': ee,
            }])],
            ignore_index=True)
        newcol = newcol + [ii]
    data['sort_index'] = newcol
    data = data.sort_values('sort_index')
    data = data.set_index('sort_index')
    name = data['name'].values
    lb = data['lower_bounds'].values
    ub = data['upper_bounds'].values
    param_bounds = list(zip(lb, ub, name))
    if debug > 2:
        print(param_bounds, '\n')
        print(N_par, '\n')
        print(name, '\n')
        print(list(map_name.keys()), '\n')
    return param_bounds, N_par, name


def read_params_results(filename):
    try:
        data = pd.read_csv(filename,
                           header=0,
                           index_col=0,
                           sep=r'\s+')
        if debug > 1:
            print(data)
    except IOError:
        func_name = inspect.currentframe().f_code.co_name
        sys.stderr.write(
            '\n ERROR in %s: cannot open file %s\nabort...\n'
            % (func_name, filename))
        sys.exit(1)
    params_list = data[0:config.NB_DIM_TOT].values
    params_listT = list(zip(*params_list))
    return params_listT


if __name__ == '__main__':
    main()
