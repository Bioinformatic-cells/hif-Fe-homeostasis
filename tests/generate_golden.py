#!/usr/bin/env python3
"""Generate / refresh golden fingerprints for HIF regression tests.

Usage (from project root):
    python tests/generate_golden.py
    python tests/generate_golden.py tests/golden/case_log_iter1.json
"""
from __future__ import print_function

import json
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

import config
import model
import simulate as sim


def _quiet():
    config.DEBUG_S = 0
    config.DEBUG_M = 0
    sim.DEBUG = 0
    model.DEBUG = 0


def generate_one(path):
    with open(path, 'r', encoding='utf-8') as fh:
        case = json.load(fh)
    _quiet()
    fingerprint = sim.evaluate_parameter_set(case['sampled_params'])
    case['expected'] = fingerprint
    with open(path, 'w', encoding='utf-8', newline='\n') as fh:
        json.dump(case, fh, indent=2, sort_keys=False)
        fh.write('\n')
    print('updated', path, 'stage=', fingerprint.get('stage'),
          'sat=', fingerprint.get('sat'))


def main(argv):
    golden_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'golden')
    if len(argv) > 1:
        paths = argv[1:]
    else:
        paths = [
            os.path.join(golden_dir, name)
            for name in sorted(os.listdir(golden_dir))
            if name.endswith('.json')
        ]
    for path in paths:
        generate_one(path)


if __name__ == '__main__':
    main(sys.argv)
