#!/usr/bin/env python3
# encoding: utf-8

import os
import sys
import tempfile
import argparse
import subprocess
from Bio import SeqIO

import recognise
from recognise.__init__ import __version__

def parse_input_args():

    parser = argparse.ArgumentParser(
        prog='recognise',
        description='recognise: identifying recombinations in genomic data',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    iGroup = parser.add_argument_group('Input options')
    iGroup.add_argument('--recipient',         '-r', help='Name of FASTA-format recipient genome',
                                                     required = True)
    iGroup.add_argument('--recombinant-list',  '-l', help='Two column list of recombinant names and FASTA-format files',
                                                     required = True)
    iGroup.add_argument('--donor',             '-d', help='Name of FASTA-format donor genome',
                                                     default = None)
    iGroup.add_argument('--version',                 action='version', version='%(prog)s {version}'.format(version=__version__))
    iGroup.add_argument('--threads',           '-c', help='Number of threads to use for parallelisation',
                                                  type=int,  default=1)
    iGroup.add_argument('--verbose',           '-v', help='Turn on debugging', action='store_true')
    iGroup.add_argument('--no-cleanup',        '-n', help='Do not cleanup intermediate files', action='store_true')

    mGroup = parser.add_argument_group('Data processing options')
    mGroup.add_argument('--method',                  help='Method to use for detecting recombinations',
                                                     choices=['3seq', 'gubbins'],
                                                     default = '3seq')

    oGroup = parser.add_argument_group('Output options')
    oGroup.add_argument('--output',            '-o', help='Output prefix',
                                                     default='recognise')

    return parser

def main():

    # import functions
    from recognise.validate import validate_input_files
    from recognise.alignment import compare_donor_recipient
    from recognise.recombination import compare_recombinant_recipient, Recombination
    from recognise.output import output_analyses

    # parse input file
    parser = parse_input_args()
    args = parser.parse_args()
    
    # set up temporary directory
    temp_dir = tempfile.TemporaryDirectory()
    sys.stderr.write('Processing data in ' + temp_dir.name + '\n')
    
    # validate input files
    recipient_id,recipient_length,recombinants = validate_input_files(args)
    
    # compare the structure of the donor and recipient
    structural_comparison = {}
    donor_maf = None
    if args.donor is not None:
        structural_comparison = compare_donor_recipient()

    # identify recombinations in recombinants
    recombinations = dict((recombinant,[]) for recombinant in recombinants.keys())
    for recombinant_name in recombinants.keys():
        recombinations[recombinant_name] = compare_recombinant_recipient(args.recipient,
                                                        recipient_id,
                                                        recombinant_name,
                                                        recombinants[recombinant_name],
                                                        args.donor,
                                                        donor_maf,
                                                        args.output,
                                                        temp_dir.name,
                                                        args.method)
                                                        
    # print summary of analysis
    check_complete = output_analyses(recombinations,recipient_length,args.output)
    # use temp_dir, and when done:
    temp_dir.cleanup()
    return check_complete

if __name__ == '__main__':
    check_val = main()
    sys.exit(check_val)