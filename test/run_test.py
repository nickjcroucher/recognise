#!/usr/bin/env python

import os
import sys
import subprocess

def clean_output(prefix):
    for suffix in []:
        os.remove('prefix' + 'suffix')

if not os.path.isdir('./data'):
    subprocess.check_output('tar xfz data.tar.gz',
                        shell=True)

subprocess.check_output('python ../recognise.py --recipient ./data/AE007317.1.fasta' + \
                        ' --recombinant-list ./data/one_recombinant.list --output ' + \
                         'one_recombinant_no_donor_test',
                        shell=True)

subprocess.check_output('python ../recognise.py --recipient ./data/AE007317.1.fasta' + \
                        ' --recombinant-list ./data/three_recombinants.list --output ' + \
                         'multiple_recombinants_no_donor_test',
                        shell=True)

subprocess.check_output('python ../recognise.py --recipient ./data/AE007317.1.fasta' + \
                        ' --recombinant-list ./data/three_recombinants.list --output ' + \
                         'one_multifasta_recombinant_no_donor_test',
                        shell=True)
