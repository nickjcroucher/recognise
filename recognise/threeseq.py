#!/usr/bin/env python3
# encoding: utf-8

import os
import sys
import re
import subprocess
from Bio import SeqIO
from Bio.Seq import Seq

def mask_alignment(alignment_file,breakpoints_3seq,recombinant_id,output_alignment_file):
    
    aln = SeqIO.parse(alignment_file, 'fasta')
    new_seqs = []

    for seq in aln:
        if seq.id == recombinant_id:
            sequence_string = str(seq.seq)
            for breakpoints in breakpoints_3seq:
                start = breakpoints[1]
                end = breakpoints[2]
                sequence_string = sequence_string[0:start] + '-'*(end-start) + sequence_string[end:len(sequence_string)]
                seq.seq = Seq(sequence_string)
        new_seqs.append(seq)
    
    SeqIO.write(new_seqs,output_alignment_file, 'fasta')
    
    return 0

def get_3seq_breakpoints(aln_file,recipient_id,donor_id,recombinant_id,length_min,threeseq_p):
    
    # wrap in try/except because 3seq reaches a threshold at which it exits with an error code
    breakpoints = []
    try:
        output_3seq = subprocess.run('yes | 3seq -triplet ' + aln_file + ' ' + \
                                                recipient_id + ' ' + \
                                                donor_id + ' ' + \
                                                recombinant_id + \
                                                ' -L' + str(length_min) + \
                                                ' -t' + str(threeseq_p) + \
                                                ' -id ' + recombinant_id,
                                                capture_output=True,
                                                shell = True)
        breakpoint_pattern = re.compile("([0-9]+)-([0-9]+) & ([0-9]+)-([0-9]+)")
        for line in output_3seq.stdout.decode().split('\n'):
            breakpoint_search = breakpoint_pattern.search(line)
            if breakpoint_search is not None:
                breakpoints.append([int(breakpoint_search.group(x)) for x in range(1,5)])
    except:
        breakpoints = []
        pass
        
    return breakpoints

def identify_recombinations_with_3seq(alignment_file,recipient_id,donor_id,recipient_mapping,tmp,prefix,recombinant,recombinant_id):

    from recognise.recombination import Recombination

    breakpoints_3seq = []
    working_alignment = os.path.join(os.path.dirname(alignment_file),'masked.'+os.path.basename(alignment_file))
    subprocess.check_output('cp ' + alignment_file + ' ' + working_alignment, shell = True)
    
    # first check if there is evidence of recombination
    continue_search = True
    
    while continue_search:
        new_breakpoints = get_3seq_breakpoints(working_alignment,recipient_id,donor_id,recombinant_id)
        if len(new_breakpoints) == 0:
            continue_search = False
        else:
            breakpoints_3seq.extend(new_breakpoints)
            mask_alignment(working_alignment,breakpoints_3seq,recombinant_id,working_alignment)
    
    # convert to recombination objects
    recombination_segments = [Recombination(start = recipient_mapping[x[1]], end = recipient_mapping[x[2]], recombinant = recombinant_id) for x in breakpoints_3seq]
    
    return recombination_segments
