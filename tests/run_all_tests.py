#!/usr/bin/env python
"""Run all regression / contract tests.

From Simulate_Py2.7_version_19_05_20:
    python tests/run_all_tests.py
"""
from __future__ import print_function

import os
import sys
import unittest

TESTS_DIR = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(TESTS_DIR)
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)


def main():
    suite = unittest.defaultTestLoader.discover(TESTS_DIR, pattern='test_*.py')
    result = unittest.TextTestRunner(verbosity=2).run(suite)
    sys.exit(0 if result.wasSuccessful() else 1)


if __name__ == '__main__':
    main()
