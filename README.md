# recognise
Precise inference of bacterial recombinations

 [![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

[![codecov](https://codecov.io/gh/nickjcroucher/recognise/branch/main/graph/badge.svg?token=Q8526GPN89)](https://codecov.io/gh/nickjcroucher/recognise)

**recognise** (**reco**mbinations in **g**e**n**om**i**c **se**quences) has been developed for the identification of recombinations in haploid genome sequences where the recipient (and donor) are known.

Recombinations can be inferred with either [Gubbins](https://github.com/nickjcroucher/gubbins), when the donor is unknown, or [3seq](https://gitlab.com/lamhm/3seq/), where the donor is known.

The input is either a two-column list of de novo assemblies (sample name and file path), or an alignment of sequences to the recipient (all in FASTA format). The recipient genome sequence is essential.

To use a list of recombinant assembly names, an example command is:
```
recognise --recipient ./data/AE007317.1.fasta --donor ./data/FQ312027.1.fasta --recombinant-list ./data/three_recombinants.list --output one_multifasta_recombinant_no_donor_test --method 3seq
```

To use an input alignment, an example command is:
```
recognise --recipient ./data/AE007317.1.fasta --donor ./data/FQ312027.1.fasta --aln-input ./data/recipient_donor_recombinant.aln --output one_multifasta_recombinant_no_donor_test  --method 3seq
```

The software can be installed through cloning the repository and creating a conda environment:

```
git clone https://github.com/nickjcroucher/recognise
conda env create -f environment.yml
cd recognise
python setup.py install
```

*This software is under active testing and development, please do not trust it for publication-level results at the moment.*
