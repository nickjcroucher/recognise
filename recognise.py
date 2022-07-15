#!/usr/bin/env python3
# encoding: utf-8

import os
import sys
import argparse

__version__ = "0.0.1"

def parse_input_args():

    parser = argparse.ArgumentParser(
        prog='recognise',
        description='recognise: identifying recombinations in genomic data',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    iGroup = parser.add_argument_group('Input options')
    iGroup.add_argument('--recipient',         '-r', help='Name of FASTA-format recipient genome', required = True)
    iGroup.add_argument('--recombinant-list',  '-l', help='List of FASTA-format recombinant genomes', required = True)
    iGroup.add_argument('--donor',             '-d', help='Name of FASTA-format donor genome')
    iGroup.add_argument('--version',                 action='version', version='%(prog)s {version}'.format(version=__version__))
    iGroup.add_argument('--threads',           '-c', help='Number of threads to use for parallelisation',
                                                  type=int,  default=1)
    iGroup.add_argument('--verbose',           '-v', help='Turn on debugging', action='store_true')
    iGroup.add_argument('--no-cleanup',        '-n', help='Do not cleanup intermediate files', action='store_true')

    mGroup = parser.add_argument_group('Data processing options')
    mGroup.add_argument('--method',                  help='Method to use for detecting recombinations',
                                                     choices=['3seq', 'gubbins'])

    oGroup = parser.add_argument_group('Output options')
    oGroup.add_argument('--output',            '-o', help='Output prefix',
                                                     default='recognise')

    return parser

def compare_donor_recipient():
    return {}

def compare_recombinant_recipient():
    return {}

def output_analyses():
    return 0

def main():
    parser = parse_input_args()
    parser.parse_args()
    # compare the structure of the donor and recipient
    structural_comparison = {}
    structural_comparison = compare_donor_recipient()
    # identify recombinations in recombinants
    recombinations = {}
    recombinations = compare_recombinant_recipient()
    # print summary of analysis
    check_complete = output_analyses()
    return check_complete

if __name__ == '__main__':
    check_val = main()
    sys.exit(check_val)
