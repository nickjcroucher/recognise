#!/usr/bin/env python3
# encoding: utf-8

import os
import sys
import tempfile
import argparse
import subprocess
from Bio import SeqIO

__version__ = "0.0.1"

class Recombination(object):
    def __init__(self,
                    start = 0,
                    end = 0,
                    p_val = 1,
                    snp_count = 0,
                    recombinant = None,
                    insertions_spanned = 0,
                    insertions_matched = 0,
                    deletions_spanned = 0,
                    deletions_matched = 0):
        self.start = start
        self.end = end
        self.p_val = p_val
        self.snp_count = snp_count
        self.recombinant = recombinant
        self.insertions_spanned = insertions_spanned
        self.insertions_matched = insertions_matched
        self.deletions_spanned = deletions_spanned
        self.deletions_matched = deletions_matched

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
    # check FASTA files
    id,length = validate_fasta(args.recipient, report = True)
    if args.donor:
        validate_fasta(args.donor)
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
    return id,length,recomb_dict

def compare_donor_recipient():
    # minigraph -cxggs R6x_GCA_000007045.1.fasta OXC141.dna > recipient_donor.gfa
    # gfatools gfa2fa serotype3_pair.gfa > serotype3_pair_merged.fa
    # gfatools bubble serotype3_pair.gfa > serotype3_pair.bed
    return {}

def identify_recombinations_with_gubbins(aln,recipient_id,recipient_mapping,tmp,prefix,recombinant,recombinant_name):
    # Get ID of recombinant - note the FASTA ID does not have to match the recombinant name
    aln_index = SeqIO.index(aln,'fasta')
    recombinant_id = None
    for id in aln_index.keys():
        if id != recipient_id:
            recombinant_id = id
    # Run Gubbins
    subprocess.check_output('run_gubbins.py --prefix ' + os.path.join(tmp,prefix) + ' --pairwise ' + aln + ' --outgroup ' + recombinant_id,
                            shell = True)
    # Process output
    recombination_list = []
    with open(os.path.join(tmp,prefix + '.recombination_predictions.gff'),'r') as rec_file:
        for line in rec_file.readlines():
            if not line.startswith('#'):
                info = line.rstrip().split()
                extra_vals = info[8].split('"')
                recombination_list.append(
                    Recombination(
                        start = recipient_mapping[int(info[3])],
                        end = recipient_mapping[int(info[4])],
                        p_val = float(extra_vals[3]),
                        snp_count = int(extra_vals[7]),
                        recombinant = recombinant_name,
                        insertions_spanned = 0,
                        insertions_matched = 0,
                        deletions_spanned = 0,
                        deletions_matched = 0
                    )
                )
    return recombination_list

def get_alignment_position_index(aln):
    aln_index = SeqIO.index(aln,'fasta')
    recipient_id = list(aln_index.keys())[0]
    recipient_sequence = aln_index[recipient_id]
    recipient_mapping = [x for x in range(0,len(recipient_sequence)+2)]
    i = 1
    j = 1
    for b in recipient_sequence:
        i = i + 1
        if b != '-':
            j = j + 1
        recipient_mapping[i] = j
    return recipient_mapping

def compare_recombinant_recipient(recipient, recipient_id, recombinant_name, recombinant, donor, prefix, tmp, method):
    # select method
    if donor is None:
        method = 'gubbins'
    # Pairwise alignment and Gubbins analysis
    if method == 'gubbins':
        # align recombinant to recipient
        raw_maf_file = os.path.join(tmp,prefix) + '.raw.maf'
        processed_maf_file = os.path.join(tmp,prefix) + '.maf'
        alignment_file = os.path.join(tmp,prefix) + '.aln'
        # Align recombinant and recipient to generate maf
        subprocess.check_output('lastz ' + recipient + ' ' + recombinant + ' --format=maf --output=' + raw_maf_file,
                                shell = True)
        # Filter maf to a single path
        subprocess.check_output('single_cov2 ' + raw_maf_file + ' > ' + processed_maf_file,
                                shell = True)
        # Convert maf to FASTA
        subprocess.check_output('maf2fasta ' + recipient + ' ' + processed_maf_file + ' fasta > ' + alignment_file,
                                shell = True)
        # Get alignment position index
        recipient_mapping = get_alignment_position_index(alignment_file)
        # Analyse output with Gubbins
        # run_gubbins.py --prefix serotype_3_mutant_comparison_og --pairwise serotype_3_mutant.unchained.aln --outgroup contig00001
        recombinations  = identify_recombinations_with_gubbins(alignment_file,
                                                                recipient_id,
                                                                recipient_mapping,
                                                                tmp,
                                                                prefix,
                                                                recombinant,
                                                                recombinant_name)
    
    # minigraph -cx lr  recipient_donor.gfa contigs.fa > recombinant.gaf
    
    return recombinations

def output_analyses(recombinations,recipient_length,prefix):
    # generate star tree and record nodes for GFF
    # output in GFF format like Gubbins for visualisation
    with open(prefix + '.recombination_predictions.gff','w') as rec_gff:
        rec_gff.write('##gff-version 3\n##sequence-region SEQUENCE 1 ' + str(recipient_length) + '\n') # need to extract recipient length
        for name in recombinations:
            for recombination in recombinations[name]:
                rec_gff.write('SEQUENCE\tGUBBINS\tCDS\t' + \
                    str(recombination.start) + '\t' + str(recombination.end) + '\t0.000\t.\t0\t' + \
                    'node="Unknown";neg_log_likelihood="' + str(recombination.p_val) + '";taxa="' + name + '";snp_count="' + str(recombination.snp_count) + '"\n');
    return 0

def main():
    parser = parse_input_args()
    args = parser.parse_args()
    # set up temporary directory
    temp_dir = tempfile.TemporaryDirectory()
    sys.stderr.write('Processing data in ' + temp_dir.name + '\n')
    # validate input files
    recipient_id,recipient_length,recombinants = validate_input_files(args)
    # compare the structure of the donor and recipient
    structural_comparison = {}
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
