#!/usr/bin/env python3
# encoding: utf-8

import os
import sys
import argparse
from Bio import SeqIO

__version__ = "0.0.1"

def parse_input_args():

    parser = argparse.ArgumentParser(
        prog='recognise',
        description='recognise: identifying recombinations in genomic data',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    iGroup = parser.add_argument_group('Input options')
    iGroup.add_argument('--recipient',         '-r', help='Name of FASTA-format recipient genome', required = True)
    iGroup.add_argument('--recombinant-list',  '-l', help='Two column list of recombinant names and FASTA-format files', required = True)
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

def validate_fasta(fn):
    if os.path.exists(fn):
        seq = SeqIO.read(fn, 'fasta')
        if '|' in seq.id:
            sys.stderr.write('Please avoid "|" in sequence IDs: ' + seq.id + ' in ' + fn + '\n')
            sys.exit(1)
    else:
        sys.stderr.write('File ' + fn + ' cannot be found\n')
        sys.exit(1)

def validate_input_files(args):
    # check FASTA files
    validate_fasta(args.recipient)
    if args.donor:
        validate_fasta(args.donor)
    # check recombinant files exist
    recomb_list = {}
    if os.path.exists(args.recombinant_list):
        with open(args.recombinant_list,'r') as recomb_file:
            for line in recomb_file:
                info = line.rstrip().split()
                if len(info) == 2:
                    validate_fasta(info[1])
                    recomb_list[info[0]] = info[1]
                else:
                    sys.stderr.write('Line does not have 2 elements: ' + line + '\n')
    else:
        sys.stderr.write('Recombinant list file ' + fn + ' cannot be found\n')
        sys.exit(1)
    return recomb_list

def compare_donor_recipient():
    # minigraph -cxggs R6x_GCA_000007045.1.fasta OXC141.dna > recipient_donor.gfa
    # gfatools gfa2fa serotype3_pair.gfa > serotype3_pair_merged.fa
    # gfatools bubble serotype3_pair.gfa > serotype3_pair.bed
    return {}

def compare_recombinant_recipient():
    # lastz R6x_GCA_000007045.1.fasta contigs.fa --format=maf --output=serotype_3_mutant.unchained.maf
    
    # single_cov2 serotype_3_mutant.unchained.maf >  serotype_3_mutant.filtered.unchained.maf
    
    # maf2fasta R6x_GCA_000007045.1.fasta serotype_3_mutant.filtered.unchained.maf fasta > serotype_3_mutant.unchained.aln
    # or
    # minigraph -cx lr  recipient_donor.gfa contigs.fa > recombinant.gaf
    
    # run_gubbins.py --prefix serotype_3_mutant_comparison_og --pairwise serotype_3_mutant.unchained.aln --outgroup contig00001
    return {}

def output_analyses():
    # output in GFF format like Gubbins for visualisation
    # generate star tree
    return 0

def main():
    parser = parse_input_args()
    args = parser.parse_args()
    # set up temporary directory
#    temp_dir = tempfile.TemporaryDirectory()
    # validate input files
    recombinant_list = validate_input_files(args)
    # compare the structure of the donor and recipient
    structural_comparison = {}
    structural_comparison = compare_donor_recipient()
    # identify recombinations in recombinants
    recombinations = {}
    recombinations = compare_recombinant_recipient()
    # print summary of analysis
    check_complete = output_analyses()
    # use temp_dir, and when done:
    
    return check_complete

if __name__ == '__main__':
    check_val = main()
    sys.exit(check_val)
