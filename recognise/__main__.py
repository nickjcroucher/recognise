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
    iGroup.add_argument('--recipient',         '-r', help='FASTA-format recipient genome',
                                                     required = True)
    iGroup.add_argument('--recombinant-list',  '-l', help='Two column list of recombinant names and FASTA-format files',
                                                     default = None)
    iGroup.add_argument('--donor',             '-d', help='FASTA-format donor genome',
                                                     default = None)
    iGroup.add_argument('--aln-input',         '-a', help='Name of FASTA-format input alignment',
                                                     default = None)
    iGroup.add_argument('--version',                 action='version', version='%(prog)s {version}'.format(version=__version__))
    iGroup.add_argument('--threads',           '-c', help='Number of threads to use for parallelisation',
                                                     type=int,
                                                     default=1)
    iGroup.add_argument('--verbose',           '-v', help='Turn on debugging', action='store_true')
    iGroup.add_argument('--no-cleanup',        '-n', help='Do not cleanup intermediate files', action='store_true')

    mGroup = parser.add_argument_group('Data processing options')
    mGroup.add_argument('--method',            '-m', help='Method to use for detecting recombinations',
                                                     choices=['3seq', 'gubbins'],
                                                     default = '3seq')
    mGroup.add_argument('--window-min',        '-w', help='Minimum window size used for Gubbins analyses',
                                                     type=int,
                                                     default=100)
    mGroup.add_argument('--window-max',        '-W', help='Maximum window size used for Gubbins analyses',
                                                     type=int,
                                                     default=10000)
    mGroup.add_argument('--snps-min',          '-s', help='Minimum number of SNPs in a recombination for Gubbins analyses',
                                                     type=int,
                                                     default=2)
    mGroup.add_argument('--length-min',        '-L', help='Minimum recombintion size for 3seq analyses',
                                                     type=int,
                                                     default=1)
    mGroup.add_argument('--gubbins-p',         '-g', help='P value threshold for identifying recombinations with Gubbins',
                                                     type=float,
                                                     default=0.05)
    mGroup.add_argument('--threeseq-p',            '-t', help='P value threshold for identifying recombinations with 3seq',
                                                     type=float,
                                                     default=0.05)
                                                     
    oGroup = parser.add_argument_group('Output options')
    oGroup.add_argument('--aln-output',        '-u', help='Store alignments in specified directory',
                                                     default=None)
    oGroup.add_argument('--output',            '-o', help='Output prefix',
                                                     default='recognise')

    return parser

def main():

    # import functions
    from recognise.validate import validate_input_files, validate_input_alignment
    from recognise.alignment import compare_donor_recipient, run_lastz_alignment
    from recognise.recombination import compare_recombinant_recipient, Recombination
    from recognise.output import output_analyses

    # parse input file
    parser = parse_input_args()
    args = parser.parse_args()
    
    # set up temporary directory
    temp_dir = tempfile.TemporaryDirectory()
    sys.stderr.write('Processing data in ' + temp_dir.name + '\n')
    
    # set up directory to store alignments if requested
    if args.aln_output is not None:
        if not os.path.isdir(args.aln_output):
            os.mkdir(args.aln_output)
    
    # validate input files
    if args.aln_input is not None:
        recipient_id,recipient_length,donor_id,recombinants,donor_in_alignment = validate_input_alignment(args)
        if not donor_in_alignment and args.method == '3seq':
            sys.stderr.write('Donor with ID ' + donor_id + ' was not found in the alignment\n')
            sys.exit(1)
    else:
        if args.recombinant_list is None:
            sys.stderr.write('Either an input alignment or list of recombinants needs to be provided\n')
            sys.exit(1)
            
        recipient_id,recipient_length,donor_id,recombinants = validate_input_files(args)
    sys.stderr.write('Completed validation of input\n')
    
    # compare the structure of the donor and recipient
    structural_comparison = {}
    donor_maf = None
    if args.donor is not None:
        if args.aln_input is None:
            donor_maf = run_lastz_alignment(args.recipient,args.donor,temp_dir.name,args.output)
        structural_comparison = compare_donor_recipient(args.recipient,args.donor)
        sys.stderr.write('Completed comparison of donor and recipient\n')
    
    # identify recombinations in recombinants
    recombinations = dict((recombinant,[]) for recombinant in recombinants.keys())
    for recombinant_name in recombinants.keys():
        recombinations[recombinant_name] = compare_recombinant_recipient(args.recipient,
                                                        recipient_id,
                                                        donor_id,
                                                        recombinant_name,
                                                        recombinants[recombinant_name],
                                                        args.donor,
                                                        donor_maf,
                                                        args.aln_input,
                                                        args.output,
                                                        args.aln_output,
                                                        temp_dir.name,
                                                        args.method,
                                                        args.window_min,
                                                        args.window_max,
                                                        args.snps_min,
                                                        args.length_min,
                                                        args.gubbins_p,
                                                        args.threeseq_p)
                                                        
    # print summary of analysis
    check_complete = output_analyses(recombinations,recipient_length,args.output)
    sys.stderr.write('Completed output file processing\n')
    # use temp_dir, and when done:
    temp_dir.cleanup()
    return check_complete

if __name__ == '__main__':
    check_val = main()
    sys.exit(check_val)
