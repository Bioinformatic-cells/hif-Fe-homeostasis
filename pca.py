#! /usr/bin/env python
# -*- coding: utf-8 -*-

"""
Module containing all the functions related to computing the pca.
"""

import numpy as np
# from numpy import array, dot, mean, std, empty, argsort

# from math import sqrt

from scipy import linalg as LA
# -----------------------------------------------------------------------------


def PCA(orig_data):
    """
   Compute the normalized PCA of a list of data with example in lines, and
   features in column
   Input: datas numpy array with features in column::orig_data

   Outputs: datas transformed to eigenspace (features in column)
   dot(evecs.T,data.T).T
       -> (evecs.T): passage matrix from eigenvector to canonical base
       -> (data.T) : dimension in row, samples in columns
       -> final .T : back to samples in rows and dimension in columns
          eigenvalues in decreasing order                               ::evals
          numpy array of eigenvector (in columns) following the
          eigenvalues order                                             ::evecs
"""
    m, n = orig_data.shape
    data = orig_data.copy()  # Copy to avoid aliasing of orig_data
    data -= np.mean(orig_data, 0)  # Center the datas
    data /= np.std(data, 0)  # Normalization of the datas with std
#    print 'norm and centered data = ',data
    R = np.cov(data, rowvar=False)  # calculate the covariance matrix
# calculate eigenvectors & eigenvalues of the covariance matrix
# use 'eigh' rather than 'eig' since R is symmetric,
# the performance gain is substantial
    evals, evecs = LA.eigh(R)
    idx = np.argsort(evals)[::-1]  # sort eigenvalue in decreasing order
    evals = evals[idx]
    evecs = evecs[:, idx]  # sort eigenvectors according to same index
    return np.dot(evecs.T, data.T).T, evals, evecs
    # np.linalg.inv(evecs) ~= evecs.T
# -----------------------------------------------------------------------------
