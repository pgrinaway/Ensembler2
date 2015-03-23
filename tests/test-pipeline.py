__author__ = 'Patrick B. Grinaway'

from Ensembler2 import search, modelling, solvation, refinement
from Ensembler2.MSMSeed import MDSys, MSMSeed
import os
import gzip
import simtk.openmm as openmm
import StringIO
def local_blast():
    """
    Test the ability to run BLAST locally.

    """
    os.environ['PDB_HOME'] = '/Users/grinawap/Documents/all_pdb/all_pdb'
    path = os.path.join('resources','4L8G_A.fasta')
    fasta_file = open(path,'r')
    fasta_text = "".join(fasta_file.readlines())
    test_msmseed = search.blast_pdb_local(fasta_text,num_hits=2)
    return test_msmseed[1]

def make_model():
    """
    Run the example through modeller
    """
    print('blasting')
    msmseed = local_blast()
    print('aligning')
    aln_seed = modelling.target_template_alignment(msmseed)
    print('modelling beginning')
    return modelling.make_model(aln_seed)

def implicit_refine_model(msmseed):
    """
    Run the example model through a 1 step refinement
    """
    print('implicit refinement beginning')
    return refinement.refine_implicitMD(msmseed,niterations=1, nsteps_per_iteration=5)

def solvate_model(msmseed):
    """
    Solvate the model, taking target_nwaters = model_nwaters+200
    """
    print('determining minimum solvation')
    nwatered_msmseed = solvation.solvate_models(msmseed)
    print('solvating to nwaters + 200')
    return solvation.solvate_models_to_target(nwatered_msmseed, target_nwaters=nwatered_msmseed.nwaters+200)

def explicit_refine_model(msmseed):
    """
    Run a very brief explicit refinement of the solvated model
    """
    print('running brief explicit refinement')
    return refinement.refine_explicitMD(msmseed)

def simulate_from_completed_models(msmseed):
    """
    This tests that the system, state, integrator are being correctly serialized.
    :param msmseed:
    :return:
    """
    serialized_system = "".join(gzip.GzipFile(fileobj=StringIO.StringIO(msmseed.explicit_refined_system)).readlines())
    serialized_state = "".join(gzip.GzipFile(fileobj=StringIO.StringIO(msmseed.explicit_refined_state)).readlines())
    serialized_integrator = "".join(gzip.GzipFile(fileobj=StringIO.StringIO(msmseed.explicit_refined_integrator)).readlines())
    system = openmm.XmlSerializer.deserialize(serialized_system)
    integrator = openmm.XmlSerializer.deserialize(serialized_integrator)
    state = openmm.XmlSerializer.deserialize(serialized_state)
    context = openmm.Context(system, integrator)
    context.setState(state)
    for i in range(5):
        print("Advancing integrator 1 step")
        integrator.step(1)
    return 0

def test_pipeline():
    msmseed = make_model()
    #print(msmseed.__dict__)
    msmeed_imp = implicit_refine_model(msmseed)
    msmeed_solv = solvate_model(msmeed_imp)
    msmseed_expl = explicit_refine_model(msmeed_solv)
    simulate_from_completed_models(msmseed_expl)
    return 0

if __name__=="__main__":
    test_pipeline()