#!/usr/bin/env python3
# encoding: utf-8

import os
import sys

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
