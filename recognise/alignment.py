#!/usr/bin/env python3
# encoding: utf-8

import os
import sys
import subprocess
from Bio import SeqIO

def run_lastz_alignment(seq1,seq2,tmp,prefix):
    seq1_fn = os.path.basename(seq1)
    seq2_fn = os.path.basename(seq2)
    raw_maf_file = os.path.join(tmp,prefix) + '.' + seq1_fn + '.' + seq2_fn + '.raw.maf'
    processed_maf_file = os.path.join(tmp,prefix) + '.' + seq1_fn + '.' + seq2_fn + '.maf'
    # Align recombinant and recipient to generate maf
    subprocess.check_output('lastz ' + seq1 + ' ' + seq2 + ' --format=maf --output=' + raw_maf_file,
                            shell = True)
    # Filter maf to a single path
    subprocess.check_output('single_cov2 ' + raw_maf_file + ' > ' + processed_maf_file,
                            shell = True)
    # Return output file name
    return processed_maf_file

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

def compare_donor_recipient(recipient,donor):
    # minigraph -cxggs R6x_GCA_000007045.1.fasta OXC141.dna > recipient_donor.gfa
    # gfatools gfa2fa serotype3_pair.gfa > serotype3_pair_merged.fa
    # gfatools bubble serotype3_pair.gfa > serotype3_pair.bed
    
    # run lastz alignment
#    donor_recipient_aln_file = run_lastz_alignment(recipient,donor)
    return {}
