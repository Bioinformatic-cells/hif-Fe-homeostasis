#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec  9 14:56:23 2019

@author: piazzai
"""

# Import the necessary packages and modules
import matplotlib.pyplot as plt
import pandas as pd
import config as config
import sys
import inspect

path = sys.argv[1] if len(sys.argv) > 1 else None


def main():
    plot_results(path)


def plot_results(path):
    data_all = []
    N_sol = int(input('How many solutions do you have?'))
    for ii in range(0, len(config.VAR_NAMES)):
        filename = path + '/results_y' + config.VAR_NAMES[ii] + '_250_' + '.txt'
        names = [config.VAR_NAMES[ii] + '__250_' + str(jj)
                 for jj in range(1, N_sol + 1)]
        try:
            data = pd.read_csv(filename,
                               usecols=['time'] + names,
                               delim_whitespace=True,
                               index_col=0,
                               engine='python')
        except IOError:
            func_name = inspect.currentframe().f_code.co_name
            print('\n ERROR in %s: cannot open file %s\nabort...'
                  % (func_name, filename), file=sys.stderr)
            sys.exit(1)
        x = data.index.values / 3600
        fig = plt.figure()
        for jj in range(0, N_sol):
            y = data[names[jj]].values
            plt.xlabel('time (h)')
            plt.ylabel('Concentration (M)')
            plt.plot(x, y)
        plt.savefig('Plot_y' + config.VAR_NAMES[ii] + '.png', dpi=100)
        plt.close(fig)
        print(data)
        data_all.append(data)
    return data_all


if __name__ == '__main__':
    main()
