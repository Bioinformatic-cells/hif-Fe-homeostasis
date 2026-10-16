#!/usr/bin/env python3
"""Regression tests for the HIF–iron homeostasis pipeline.

These tests freeze the numerical fingerprint of known parameter sets.
If a code change alters the biology/numerics, the test fails.

Refresh goldens only after an intentional, reviewed change:
    python tests/generate_golden.py
"""
from __future__ import print_function

import json
import os
import sys
import unittest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

import numpy as np

import config
import model
import simulate as sim

GOLDEN_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'golden')
RTOL = 1e-6
ATOL = 1e-12


def _quiet():
    config.DEBUG_S = 0
    config.DEBUG_M = 0
    sim.DEBUG = 0
    model.DEBUG = 0


def _load_cases():
    cases = []
    for name in sorted(os.listdir(GOLDEN_DIR)):
        if not name.endswith('.json'):
            continue
        path = os.path.join(GOLDEN_DIR, name)
        with open(path, 'r', encoding='utf-8') as fh:
            case = json.load(fh)
        case['_path'] = path
        cases.append(case)
    return cases


def _assert_close(test, got, expected, label):
    if isinstance(expected, (list, tuple)):
        test.assertEqual(len(got), len(expected), label + ' length')
        for i, (g, e) in enumerate(zip(got, expected)):
            _assert_close(test, g, e, '%s[%d]' % (label, i))
        return
    if isinstance(expected, bool):
        test.assertEqual(bool(got), expected, label)
        return
    if isinstance(expected, str):
        test.assertEqual(got, expected, label)
        return
    if expected is None:
        test.assertIsNone(got, label)
        return
    test.assertTrue(
        np.isclose(float(got), float(expected), rtol=RTOL, atol=ATOL),
        '%s: got %s expected %s' % (label, got, expected))


class TestHIFRegression(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        _quiet()
        cls.cases = _load_cases()

    def test_golden_cases_present(self):
        self.assertGreater(len(self.cases), 0, 'no golden JSON files found')

    def test_parameter_sets_match_golden_fingerprints(self):
        for case in self.cases:
            expected = case.get('expected')
            self.assertIsNotNone(
                expected,
                'missing expected fingerprint in %s; run '
                'python tests/generate_golden.py' % case['_path'])
            got = sim.evaluate_parameter_set(case['sampled_params'])
            self.assertEqual(got['stage'], expected['stage'], case['name'])
            self.assertEqual(got['sat'], expected['sat'], case['name'])
            self.assertEqual(got['err'], expected['err'], case['name'])
            for key in ('algebraic_init_vars', 'init_vars', 'end_states',
                        'relative_sols', 'Fe_tot', 'params'):
                if key in expected and expected[key]:
                    _assert_close(self, got.get(key, []), expected[key],
                                  '%s.%s' % (case['name'], key))


if __name__ == '__main__':
    unittest.main()
