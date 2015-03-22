__author__ = 'Patrick B. Grinaway'

from Ensembler2 import MSMSeed


def _PIR_alignment(target_sequence, target_id, template_sequence, template_id):
    """
    Produce a pairwise alignment of the target sequences and template sequence in the format
    that MODELLER likes.

    Arguments
    ---------
    target_sequence : String
        The sequence of the target protein in FASTA format
    target_id : String
        The name of the target
    template_sequence : String
        The sequence of the target protein in FASTA format
    template_id : String
        The name of the template

    Returns
    -------
    contents : String
        The result of the sequence alignment in PIR format (for MODELLER)
    """
    import Bio.SubsMat
    import Bio.pairwise2
    from Bio.SubsMat import MatrixInfo as matlist
    import Bio.SeqIO

    matrix = matlist.gonnet
    gap_open = -10
    gap_extend = -0.5
    aln = Bio.pairwise2.align.globalds(target_sequence, template_sequence, matrix, gap_open, gap_extend)

    #put together PIR file
    contents = "Target-template alignment by clustal omega\n"
    contents += ">P1;%s\n" % target_id
    contents += "sequence:%s:FIRST:@:LAST :@:::-1.00:-1.00\n" % target_id
    contents += aln[0][0] + '*\n'
    contents += ">P1;%s\n" % template_id
    contents += "structureX:%s:FIRST:@:LAST : :undefined:undefined:-1.00:-1.00\n" % template_id
    contents += aln[0][1] + '*\n'
    return contents

def target_template_alignment(msmseed):
    """
    Use Biopython pairwise2 to generate a sequence alignment in PIR format so that MODELLER can be used.
    Puts alignment in MSMSeed. The functionality for this was separated into another internal function,
    _PIR_alignment().

    Parameters
    ----------

    msmseed : MSMSeed
        object containing the sequences to be aligned

    Returns
    -------
        msmseed : MSMSeed
            Object containing the alignment between the input sequences accessible from the alignment property


    """
    aln = _PIR_alignment(msmseed.target_sequence, msmseed.target_id, msmseed.template_sequence, msmseed.template_id)
    msmseed.alignment = str(aln)
    return msmseed


def align_template_to_reference(msmseed, ref_msmseed):
    """
    Performs a structural alignment between the reference and template
    to retrieve an RMSD

    Arguments
    ---------
    msmseed : MSMSeed object
        An MSMSeed object containing the template to align
    ref_msmseed : MSMSeed object
        An MSMSeed object containing the reference structure

    Returns
    -------
    msmseed : MSMSeed object
        The template MSMSeed object with an RMSD to template attribute
    """
    import modeller
    import tempfile
    import shutil
    import copy
    import os
    temp_dir = tempfile.mkdtemp()
    try:
        os.chdir(temp_dir)
        alignment_file = open('aln_tmp.pir','w')
        aln = _PIR_alignment(ref_msmseed.template_sequence, ref_msmseed.template_id, msmseed.template_sequence, msmseed.template_id)
        alignment_file.writelines(aln)
        alignment_file.close()
        template_file = open(msmseed.template_id + '.pdb','w')
        template_pdb = msmseed.template_structure
        template_pdb.writeFile(template_pdb.topology, template_pdb.positions, template_file)
        template_file.close()
        ref_pdb = ref_msmseed.template_structure
        ref_file = open(ref_msmseed.template_id + '.pdb', 'w')
        ref_pdb.writeFile(ref_pdb.topology, ref_pdb.positions, ref_file)
        ref_file.close()
        modeller.log.none()
        env = modeller.environ()
        env.io.atom_files_directory = temp_dir
        aln = modeller.alignment(env, file='aln_tmp.pir', align_codes=(ref_msmseed.template_id, msmseed.template_id))
        mdl  = modeller.model(env, file=ref_msmseed.template_id + '.pdb')
        mdl2 = modeller.model(env, file=msmseed.template_id+'.pdb')
        atmsel = modeller.selection(mdl).only_atom_types('CA')
        r = atmsel.superpose(mdl2, aln)
        msmseed.rmsd_to_reference = copy.deepcopy(r.rms)
    except Exception as e:
        msmseed.error_message = e.message
    finally:
        shutil.rmtree(temp_dir)
    return msmseed

def target_template_alignment(msmseed):
    """
    Use Biopython pairwise2 to generate a sequence alignment in PIR format so that MODELLER can be used.
    Puts alignment in MSMSeed. The functionality for this was separated into another internal function,
    _PIR_alignment().

    Parameters
    ----------

    msmseed : MSMSeed
        object containing the sequences to be aligned

    Returns
    -------
        msmseed : MSMSeed
            Object containing the alignment between the input sequences accessible from the alignment property


    """
    aln = _PIR_alignment(msmseed.target_sequence, msmseed.target_id, msmseed.template_sequence, msmseed.template_id)
    msmseed.alignment = str(aln)
    return msmseed

def make_model(msmseed):
    """
    Use MODELLER from the Sali lab to create a model between the target and template specified in the input

    Parameters
    ----------
    msmseed : MSMSeed
        object containing the alignment between target and template and template structure

    Returns
    -------
    msmseed : MSMSeed
        object containing the homology model built from the input alignment and template structure
    """

    import tempfile
    import os
    import modeller
    import modeller.automodel
    import shutil
    import simtk.openmm.app as app
    #if the target and template are the same, modeller dies.
    if msmseed.template_id == msmseed.target_id:
        msmseed.target_model = msmseed.template_structure
        return msmseed
    #first, we need to make a temp directory where we can put the files MODELLER needs
    temp_dir = tempfile.mkdtemp()
    try:
        os.chdir(temp_dir)
        alignment_file = open('aln_tmp.pir','w')
        alignment_file.writelines(msmseed.alignment)
        alignment_file.close()
        template_file = open(msmseed.template_id + '.pdb','w')
        template_pdb = msmseed.template_structure
        template_pdb.writeFile(template_pdb.topology, template_pdb.positions, template_file)
        template_file.close()
        modeller.log.none()
        env = modeller.environ()
        env.io.atom_files_directory = temp_dir
        a = modeller.automodel.allhmodel(env,
                                         # file with template codes and target sequence
                                         alnfile  = 'aln_tmp.pir',
                                         # PDB codes of the template
                                         knowns   = msmseed.template_id,
                                         # code of the target
                                         sequence = msmseed.target_id)
        a.make()
        tmp_model_pdbfilename = a.outputs[0]['name']
        target_model = modeller.model(env, file=tmp_model_pdbfilename)
        msmseed.sequence_similarity = target_model.seq_id
        msmseed.target_model = app.PDBFile(tmp_model_pdbfilename)
        msmseed.target_restraints = open('%s.rsr' % msmseed.target_id, 'r').readlines()
    except:
        msmseed.error_state = -1

    finally:
        shutil.rmtree(temp_dir)
    return msmseed
