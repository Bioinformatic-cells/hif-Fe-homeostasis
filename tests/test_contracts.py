#!/usr/bin/env python3
"""Lightweight contract tests that do not run the ODE solver."""
from __future__ import print_function

import os
import sys
import unittest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

import config


class TestHIFConfigContracts(unittest.TestCase):
    def test_param_names_match_dimensions(self):
        self.assertEqual(len(config.PARAM_NAMES), config.NB_DIM_TOT)
        self.assertEqual(len(config.map_PARAM_NAMES), config.NB_SAMPLED_DIM)
        self.assertEqual(
            config.NB_SAMPLED_DIM + len(config.PARAM_ALGEB_DEFS),
            config.NB_DIM_TOT)

    def test_var_names_match_dimensions(self):
        self.assertEqual(len(config.VAR_NAMES), config.NB_VARS)
        self.assertEqual(len(config.map_VAR_NAMES), config.NB_VARS)
        self.assertEqual(len(config.VARS_BOUNDS), config.NB_VARS)

    def test_constraint_box_covers_sampled_params(self):
        self.assertEqual(len(config.PARAMS_CONSTRAINT_BOX), config.NB_SAMPLED_DIM)

    def test_inputs_consistent(self):
        self.assertEqual(len(config.MULTIPLE_INPUTS), config.NB_INPUTS)
        self.assertGreaterEqual(config.NB_INPUTS, 1)

    def test_map_indices_are_1_based_unique(self):
        values = sorted(config.map_PARAM_NAMES.values())
        self.assertEqual(values, list(range(1, config.NB_SAMPLED_DIM + 1)))
        var_values = sorted(config.map_VAR_NAMES.values())
        self.assertEqual(var_values, list(range(1, config.NB_VARS + 1)))


class TestBoundFiles(unittest.TestCase):
    def test_hif_bound_files_exist_and_parse(self):
        import read_file
        vfile = os.path.join(ROOT, 'b-vars-hif.txt')
        pfile = os.path.join(ROOT, 'b-par-hif.txt')
        self.assertTrue(os.path.isfile(vfile))
        self.assertTrue(os.path.isfile(pfile))
        var_bounds, n_vars, var_names = read_file.read_file(
            vfile, 0, config.map_VAR_NAMES)
        par_bounds, n_pars, par_names = read_file.read_file(
            pfile, 0, config.map_PARAM_NAMES)
        self.assertEqual(n_vars, config.NB_VARS)
        self.assertEqual(n_pars, config.NB_SAMPLED_DIM)
        self.assertEqual(len(var_bounds), config.NB_VARS)
        self.assertEqual(len(par_bounds), config.NB_SAMPLED_DIM)


if __name__ == '__main__':
    unittest.main()
