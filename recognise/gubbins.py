#!/usr/bin/env python3
# encoding: utf-8

import os
import sys
import subprocess
from Bio import SeqIO

def identify_recombinations_with_gubbins(aln,recipient_id,recipient_mapping,tmp,prefix,recombinant,
                                            recombinant_name,window_min,window_max,snps_min,gubbins_p):

    # Avoid circular import
    from recognise.recombination import Recombination

    # Get ID of recombinant - note the FASTA ID does not have to match the recombinant name
    aln_index = SeqIO.index(aln,'fasta')
    recombinant_id = None
    for id in aln_index.keys():
        if id != recipient_id:
            recombinant_id = id
    # Run Gubbins
    subprocess.check_output('run_gubbins.py --prefix ' + os.path.join(tmp,prefix) + \
                            ' --pairwise ' + \
                            ' --min-window-size ' + str(window_min) + \
                            ' --max-window-size ' + str(window_max) + \
                            ' --p-value ' + str(gubbins_p) + \
                            ' --min-snps ' + str(snps_min) + \
                            ' --outgroup ' + recombinant_id + \
                            ' ' + aln + \
                            ' &> /dev/null',
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
