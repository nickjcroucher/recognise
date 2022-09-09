#!/usr/bin/env python

import os
import sys
import unittest
import subprocess

def clean_output(prefix):
    for suffix in ['.recombination_predictions.gff']:
        os.remove('prefix' + 'suffix')

if not os.path.isdir('./data'):
    subprocess.check_output('tar xfz data.tar.gz',
                        shell=True)

class RecogniseTests(unittest.TestCase):

    def test_one_contiguous_recombinant_no_donor(self):
        subprocess.check_output('recognise --recipient ./data/AE007317.1.fasta' + \
                                ' --recombinant-list ./data/one_recombinant.list --output ' + \
                                 'one_recombinant_no_donor_test',
                                shell=True)
                                
    def test_three_contiguous_recombinants_no_donor(self):
        subprocess.check_output('recognise --recipient ./data/AE007317.1.fasta' + \
                                ' --recombinant-list ./data/three_recombinants.list --output ' + \
                                 'multiple_recombinants_no_donor_test',
                                shell=True)
                            
    def test_one_noncontiguous_recombinant_no_donor(self):
        subprocess.check_output('recognise --recipient ./data/AE007317.1.fasta' + \
                                ' --recombinant-list ./data/three_recombinants.list --output ' + \
                                 'one_multifasta_recombinant_no_donor_test',
                                shell=True)

    def test_one_contiguous_recombinant_no_donor_store_alignment(self):
        subprocess.check_output('recognise --recipient ./data/AE007317.1.fasta' + \
                                ' --aln-output alignments ' + \
                                ' --recombinant-list ./data/one_recombinant.list --output ' + \
                                 'one_recombinant_no_donor_test',
                                shell=True)

    def test_one_noncontiguous_recombinant_with_donor(self):
        subprocess.check_output('recognise --recipient ./data/AE007317.1.fasta' + \
                                ' --donor ./data/FQ312027.1.fasta ' + \
                                ' --recombinant-list ./data/three_recombinants.list --output ' + \
                                 'one_multifasta_recombinant_no_donor_test',
                                shell=True)

    def test_one_contiguous_recombinant(self):
        subprocess.check_output('recognise --recipient ./data/AE007317.1.fasta' + \
                                ' --recombinant-list ./data/one_recombinant.list --output ' + \
                                ' --donor ./data/FQ312027.1.fasta --method 3seq ' + \
                                 'one_recombinant_test',
                                shell=True)

    def test_three_contiguous_recombinants(self):
        subprocess.check_output('recognise --recipient ./data/AE007317.1.fasta' + \
                                ' --recombinant-list ./data/three_recombinants.list --output ' + \
                                ' --donor ./data/FQ312027.1.fasta --method 3seq ' + \
                                 'multiple_recombinants_test',
                                shell=True)

    def test_one_noncontiguous_recombinant(self):
        subprocess.check_output('recognise --recipient ./data/AE007317.1.fasta' + \
                                ' --recombinant-list ./data/three_recombinants.list --output ' + \
                                ' --donor ./data/FQ312027.1.fasta --method 3seq ' + \
                                 'one_multifasta_recombinant_test',
                                shell=True)
