#!/usr/bin/env python

import argparse
from pyfaidx import Faidx

parser = argparse.ArgumentParser()
parser.add_argument('-g', '--genome', type=str, required=True,
                    help='Path to genome sequence')
args = parser.parse_args()

Faidx(args.genome)
