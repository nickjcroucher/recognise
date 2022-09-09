#!/usr/bin/env python3
# encoding: utf-8

import os
import sys
import subprocess
from Bio import SeqIO

from recognise.alignment import get_alignment_position_index, run_lastz_alignment
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
        
def compare_recombinant_recipient(recipient, recipient_id, donor_id, recombinant_name, recombinant, donor, donor_maf, prefix, aln_dir, tmp, method):
    
    # select method
    if (donor is None or donor_maf is None) and method == '3seq':
        sys.stderr.write('Unable to align donor and recipient - skipping 3seq analysis\n')
        method = 'gubbins'
    
    ######################
    # Pairwise alignment #
    ######################
    
    # run pairwise alignment of recipient and recombinant with lastz
    recombinant_maf_file = run_lastz_alignment(recipient,recombinant,tmp,prefix)
    
    ####################
    # Gubbins analysis #
    ####################

    # Convert pairwise maf to FASTA
    alignment_file = os.path.join(tmp,prefix) + os.path.basename(recombinant) + '.aln'
    subprocess.check_output('maf2fasta ' + recipient + ' ' + recombinant_maf_file + ' fasta > ' + alignment_file,
                            shell = True)
    # Get alignment position index
    recipient_mapping = get_alignment_position_index(alignment_file,recipient_id)
    # Analyse output with Gubbins
    # run_gubbins.py --prefix serotype_3_mutant_comparison_og --pairwise serotype_3_mutant.unchained.aln --outgroup contig00001
    mosaic_recombinations  = identify_recombinations_with_gubbins(alignment_file,
                                                            recipient_id,
                                                            recipient_mapping,
                                                            tmp,
                                                            prefix,
                                                            recombinant,
                                                            recombinant_name)
    
    #################
    # 3seq analysis #
    #################

    if method == '3seq':
        # combine alignments of donor and recombinant to recipient
        subprocess.check_output('multiz ' + donor_maf + ' ' + recombinant_maf_file + ' 1 > ' + recombinant_name + '.multiz',
                                shell = True)
        alignment_file_3seq = os.path.join(tmp,prefix) + '.' + recombinant_name + '.3seq.aln'
        subprocess.check_output('maf2fasta ' + recipient + ' ' + recombinant_name + '.multiz' + ' fasta > ' + alignment_file_3seq,
                                shell = True)
#        alignment_file_3seq = './alignments/one_recombinant_no_donor_test.3seq.aln' # debug
        recipient_mapping_3seq = get_alignment_position_index(alignment_file_3seq,recipient_id)
        # iteratively identify breakpoints in the triplet
        mosaic_recombinations = identify_recombinations_with_3seq(alignment_file_3seq,
                                                                    recipient_id,
                                                                    donor_id,
                                                                    recipient_mapping_3seq,
                                                                    tmp,
                                                                    prefix,
                                                                    recombinant,
                                                                    recombinant_name
        )
        
        
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
