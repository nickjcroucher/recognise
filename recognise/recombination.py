#!/usr/bin/env python3
# encoding: utf-8

import os
import sys
import subprocess
from Bio import SeqIO

from recognise.alignment import get_alignment_position_index, run_lastz_alignment, extract_subalignment
from recognise.gubbins import identify_recombinations_with_gubbins
from recognise.threeseq import identify_recombinations_with_3seq

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
        
def compare_recombinant_recipient(recipient, recipient_id, donor_id, recombinant_name, recombinant, donor,
                                    donor_maf, input_aln, prefix, aln_dir, tmp, method, window_min,
                                    window_max, snps_min, length_min, gubbins_p, threeseq_p):
    
    ######################
    # Pairwise alignment #
    ######################
    
    if input_aln is None:
        # run pairwise alignment of recipient and recombinant with lastz
        recombinant_maf_file = run_lastz_alignment(recipient,recombinant,tmp,prefix)
    
    ####################
    # Gubbins analysis #
    ####################

    # Convert pairwise maf to FASTA
    alignment_file = os.path.join(tmp,prefix) + os.path.basename(recombinant) + '.aln'
    if input_aln is None:
        subprocess.check_output('maf2fasta ' + recipient + ' ' + recombinant_maf_file + ' fasta > ' + alignment_file,
                            shell = True)
    else:
        extract_subalignment(input_aln,[recipient_id,recombinant_name],alignment_file)
    # Get alignment position index
    recipient_mapping = get_alignment_position_index(alignment_file,recipient_id)
    sys.stderr.write('Completed pairwise alignment of recombinant and recipient\n')
    # Analyse output with Gubbins
    # run_gubbins.py --prefix serotype_3_mutant_comparison_og --pairwise serotype_3_mutant.unchained.aln --outgroup contig00001
    mosaic_recombinations  = identify_recombinations_with_gubbins(alignment_file,
                                                            recipient_id,
                                                            recipient_mapping,
                                                            tmp,
                                                            prefix,
                                                            recombinant,
                                                            recombinant_name,
                                                            window_min,
                                                            window_max,
                                                            snps_min,
                                                            gubbins_p)
    sys.stderr.write('Completed analysis with Gubbins\n')
    
    #################
    # 3seq analysis #
    #################

    if method == '3seq':
        alignment_file_3seq = os.path.join(tmp,prefix) + '.' + recombinant_name + '.3seq.aln'
        if input_aln is None:
            # combine alignments of donor and recombinant to recipient
            subprocess.check_output('multiz ' + donor_maf + ' ' + recombinant_maf_file + ' 1 > ' + recombinant_name + '.multiz',
                                    shell = True)
            subprocess.check_output('maf2fasta ' + recipient + ' ' + recombinant_name + '.multiz' + ' fasta > ' + alignment_file_3seq,
                                    shell = True)
        else:
            extract_subalignment(input_aln,[recipient_id,donor_id,recombinant_name],alignment_file_3seq)
#        alignment_file_3seq = './alignments/one_recombinant_no_donor_test.3seq.aln' # debug
        recipient_mapping_3seq = get_alignment_position_index(alignment_file_3seq,recipient_id)
        sys.stderr.write('Completed alignment of recombinant, recipient and donor\n')
        # iteratively identify breakpoints in the triplet
        mosaic_recombinations = identify_recombinations_with_3seq(alignment_file_3seq,
                                                                    recipient_id,
                                                                    donor_id,
                                                                    recipient_mapping_3seq,
                                                                    tmp,
                                                                    prefix,
                                                                    recombinant,
                                                                    recombinant_name,
                                                                    length_min,
                                                                    threeseq_p)
        sys.stderr.write('Completed analysis with 3seq\n')
    # minigraph -cx lr  recipient_donor.gfa contigs.fa > recombinant.gaf
    
    ####################
    # Store alignments #
    ####################
    
    if aln_dir is not None:
        subprocess.check_output('cp ' + alignment_file + ' ' + aln_dir,
                                shell = True)
        if method == '3seq':
            subprocess.check_output('cp ' + alignment_file_3seq + ' ' + aln_dir,
                            shell = True)
    
    return mosaic_recombinations
