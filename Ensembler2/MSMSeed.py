__author__ = 'Patrick B. Grinaway'



class MSMSeed(object):
    """
    Container class for information used in creating starting structures with MSMseeder

    Parameters
    ----------

    target_sequence : String
        FASTA format string with the target amino acid sequence and name of target after >
    template_sequence : String
        FASTA format string with the template amino acid sequence and name of template after >
    template_structure : String
        contains the string representation of the template pdb structure
    blast_eval : float, optional
        contains the e-value obtained from blast when aligning the template and target. Default 0


    Attributes
    ----------
    blast_eval : float
        contains the e-value obtained via blast search.
    rmsd_to_reference : float
        contains the RMSD to a reference structure
    target_sequence : String
        containing the sequence of this model's target in fasta
    target_ID : String
        containing the name of this model's target
    template_ID : String
        containing the name of this model's template
    template_sequence : String
        containing the sequence of this model's template in fasta
    template_structure : simtk.openmm.app.PDBFile
        object containing the structure of this model's template
    alignment : String
        containing the PIR format (MODELLER style) alignment of the target and template
    target_model : simtk.openmm.app.PDBFile
        object containing the output of MODELLER
    target_restraints: String
        containing .rsr file output of MODELLER
    sequence_similarity : float
        contains the sequence similarity between the target and template
    implicit_refined_model : MDSys
        object containing an openmm topology and positions of an implicitly-refined model
    error_message : String
        containing a brief explanation of the error condition on this model, if any
    error_state : int
        the stage at which the error occurred
    nwaters : int
        the number of waters used to solvate the model
    solvated_model : MDSys
        object containing an openmm topology and positions of the solvated model
    explicit_refined_pdb : String
        gzipped PDB file contents of explicitly-refined model
    explicit_refined_state : String
        gzipped xml serialized form of explicit-refined model state
    explicit_refined_system : String
        gzipped xml serialized form of explicit-refined model system
    explicit_refined_integrator : String
        gzipped xml serialized form of explicit-refined model integrator





    """

    def __init__(self, target_sequence, template_sequence, template_structure, blast_eval=0):
        import StringIO
        import Bio.SeqIO

        target_seq_record = list()
        template_seq_record = list()

        target_sequence_stringio = StringIO.StringIO(target_sequence)
        template_sequence_stringio = StringIO.StringIO(template_sequence)

        target_sequence_stringio.seek(0)
        target_sequence_stringio.seek(0)

        target_sequence_parser = Bio.SeqIO.parse(target_sequence_stringio,'fasta')
        template_sequence_parser = Bio.SeqIO.parse(template_sequence_stringio,'fasta')

        for record in target_sequence_parser:
            self._target_sequence= record.seq
            self._target_id = record.id
        for record in template_sequence_parser:
            self._template_sequence = record.seq
            self._template_id = record.id

        self._template_structure = template_structure
        self._error_state = 0
        self._nwaters = 0
        self._blast_eval = blast_eval


    @property
    def blast_eval(self):
        return self._blast_eval

    @property
    def rmsd_to_reference(self):
        return self._rmsd_to_reference
    @rmsd_to_reference.setter
    def rmsd_to_reference(self, rmsd):
        self._rmsd_to_reference = rmsd

    @property
    def target_sequence(self):
        return self._target_sequence
    @property
    def template_sequence(self):
        return self._template_sequence
    @property
    def template_structure(self):
        return self._template_structure

    @property
    def alignment(self):
        return self._alignment
    @alignment.setter
    def alignment(self, new_alignment):
        self._alignment = new_alignment

    @property
    def target_model(self):
        return self._target_model
    @target_model.setter
    def target_model(self, new_target_model):
        self._target_model = new_target_model

    @property
    def sequence_similarity(self):
        return self._sequence_similarity
    @sequence_similarity.setter
    def sequence_similarity(self, similarity):
        self._sequence_similarity = similarity

    @property
    def target_restraints(self):
        return self._target_restraints
    @target_restraints.setter
    def target_restraints(self,new_target_restraints):
        self._target_restraints = new_target_restraints

    @property
    def target_id(self):
        return self._target_id
    @property
    def template_id(self):
        return self._template_id

    @property
    def implicit_refined_model(self):
        return self._implicit_refined_model
    @implicit_refined_model.setter
    def implicit_refined_model(self, new_model):
        self._implicit_refined_model = new_model

    @property
    def error_state(self):
        return self._error_state
    @error_state.setter
    def error_state(self, error):
        self._error_state = error

    @property
    def error_message(self):
        return self._error_message
    @error_message.setter
    def error_message(self, message):
        self._error_message = message

    @property
    def nwaters(self):
        return self._nwaters
    @nwaters.setter
    def nwaters(self,num_waters):
        self._nwaters = num_waters


    @property
    def solvated_model(self):
        return self._solvated_model
    @solvated_model.setter
    def solvated_model(self, new_solvated_model):
        self._solvated_model = new_solvated_model

    @property
    def explicit_refined_pdb(self):
        return self._explicit_refined_pdb
    @explicit_refined_pdb.setter
    def explicit_refined_pdb(self, new_pdb):
        self._explicit_refined_pdb = new_pdb

    @property
    def explicit_refined_state(self):
        return self._explicit_refined_state
    @explicit_refined_state.setter
    def explicit_refined_state(self,new_state):
        self._explicit_refined_state = new_state

    @property
    def explicit_refined_integrator(self):
        return self._explicit_refined_integrator
    @explicit_refined_integrator.setter
    def explicit_refined_integrator(self, new_integrator):
        self._explicit_refined_integrator = new_integrator

    @property
    def explicit_refined_system(self):
        return self._explicit_refined_system
    @explicit_refined_system.setter
    def explicit_refined_system(self, new_system):
        self._explicit_refined_system = new_system



class MDSys(object):
    def __init__(self, topology, positions):
        self._topology = topology
        self._positions = positions
    @property
    def topology(self):
        return self._topology
    @property
    def positions(self):
        return self._positions

