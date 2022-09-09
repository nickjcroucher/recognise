#!/usr/bin/env python3
# encoding: utf-8

import os
import sys
from Bio import SeqIO

def validate_fasta(fn, report = False):
    if os.path.exists(fn):
        seq = SeqIO.read(fn, 'fasta')
        if '|' in seq.id:
            sys.stderr.write('Please avoid "|" in sequence IDs: ' + seq.id + ' in ' + fn + '\n')
            sys.exit(1)
        if report:
            return seq.id,len(seq.seq)
    else:
        sys.stderr.write('File ' + fn + ' cannot be found\n')
        sys.exit(1)

def validate_input_files(args):
    # check FASTA file for recipient
    id,length = validate_fasta(args.recipient, report = True)
    # check FASTA file for donor
    donor_id = None
    if args.donor:
        donor_id,length = validate_fasta(args.donor, report = True)
    # check recombinant files exist
    recomb_dict = {}
    if os.path.exists(args.recombinant_list):
        with open(args.recombinant_list,'r') as recomb_file:
            for line in recomb_file:
                info = line.rstrip().split()
                if len(info) == 2:
                    validate_fasta(info[1])
                    recomb_dict[info[0]] = info[1]
                else:
                    sys.stderr.write('Line does not have 2 elements: ' + line + '\n')
    else:
        sys.stderr.write('Recombinant list file ' + fn + ' cannot be found\n')
        sys.exit(1)
    return id,length,donor_id,recomb_dict
