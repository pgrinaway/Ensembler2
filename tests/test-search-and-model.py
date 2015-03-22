__author__ = 'Patrick B. Grinaway'


import os
from Ensembler2 import MSMSeed, search, modelling

def local_blast():
    """
    Test the ability to run BLAST locally.

    """
    path = os.path.join('tests','resources','4L8G_A.fasta')
    fasta_file = open(path,'r')
    fasta_text = "".join(fasta_file.readlines())
    test_msmseed = search.blast_pdb_local(fasta_text,num_hits=1)
    return test_msmseed

def make_model():
    msmseed = local_blast()
    return modelling.make_model(msmseed)

def make_structural_alignment():


