__author__ = 'Patrick B. Grinaway'
import pyspark
import os
from Ensembler2 import modelling, MSMSeed, refinement, search, solvation
import simtk.openmm.app as app

class SparkDriver(object):
    """
    Driver class for the Spark version of Ensembler2. Provides a higher-
    level interface for interacting with Ensembler2.
    """
    _starting_templates = None
    _rdd_starting_templates = None
    _modeled_seeds = None
    _filtered_models = None
    _failed_seeds = None
    _implicit_refined_seeds = None
    _failure_data = None
    _solvated_models = None

    def __init__(self, target_squence, spark_master, models_directory, template_model_directory, constrained_memory=True, executor_memory='5g', auto_initialize=True, platform = 'CUDA', n_per_slice=10):
        self._target_sequence = target_squence
        spark_conf = pyspark.SparkConf()
        spark_conf.setMaster(spark_master)
        spark_conf.set("spark.executor.memory", executor_memory)
        self._template_model_directory = os.path.abspath(template_model_directory)
        self._spark_context = pyspark.SparkContext(conf=spark_conf)
        self._models_directory = os.path.abspath(models_directory)
        self._model_rdds = dict()
        self._platform = platform
        self._n_per_slice = n_per_slice
        self._constrained_memory = constrained_memory
        self._error_data = list()
        if auto_initialize:
            self.initialize_with_blast()



    def initialize_with_blast(self, local=True, max_hits=1000):
        """
        Initializes the SparkDriver class with results from a BLAST query, up to max_hits
        using local blast if specified.
        """
        if local:
            self._starting_templates = search.blast_pdb_local(self._target_sequence, num_hits=max_hits)
        else:
            self._starting_templates = search.blast_pdb(self._target_sequence, num_hits=max_hits)
        self._rdd_starting_templates = self._spark_context.parallelize(self._starting_templates)


    def make_models(self):
        """
        Run modeller on identified homologues, collecting failure data
        """

        self._modeled_seeds = self._rdd_starting_templates.map(modelling.target_template_alignment).map(modelling.make_model)
        self._error_data.extend(self._modeled_seeds.filter(lambda seed: seed.error_state < 0).map(self.get_error_metadata).collect()) #TODO: do something to collect failures
        self._modeled_seeds = self._modeled_seeds.filter(lambda seed: seed.error_state == 0)

    def refine_implicit(self):
        """
        Perform implicit-solvent refinement of the models
        :return:
        """
        def implicit_refine(x):
            result = refinement.refine_implicitMD(x, openmm_platform='CUDA', niterations=500, nsteps_per_iteration=100)
            return result
        self._implicit_refined_seeds = self._modeled_seeds.map(implicit_refine).persist(pyspark.StorageLevel.MEMORY_AND_DISK)
        self._error_data.extend(self._implicit_refined_seeds.filter(lambda seed: seed.error_state < 0).map(self.get_error_metadata).collect())
        self._implicit_refined_seeds = self._implicit_refined_seeds.filter(lambda seed: seed.error_state == 0).persist(pyspark.StorageLevel.MEMORY_AND_DISK)

    def solvate_models(self):
        """
        Determine the number of waters needed to minimally solvate the model
        :return:
        """
        def get_nwaters(x):
            return x.nwaters
        nwaters_array = self._implicit_refined_seeds.map(solvation.solvate_models).map(get_nwaters).collect()
        nwaters_target = solvation.calculate_nwaters(nwaters_array)

        def solvate_target(x):
            result = solvation.solvate_models_to_target(x, nwaters_target)
            return result
        self._solvated_models = self._implicit_refined_seeds.map(solvate_target)
        self._error_data.extend(self._solvated_models.filter(lambda seed: seed.error_state < 0).map(self.get_error_metadata).collect())
        self._solvated_models = self._solvated_models.filter(lambda seed: seed.error_state == 0).persist(pyspark.StorageLevel.MEMORY_AND_DISK)

    def explicit_refine_models(self):
        """
        Perform explicit refinement of models
        :return:
        """
        def explicit_refine(x):
            result = refinement.refine_explicitMD(x, openmm_platform='CUDA', niterations=500, nsteps_per_iteration=100)
            return result
        self._explicit_refined_models = self._solvated_models.map(explicit_refine).persist(pyspark.StorageLevel.MEMORY_AND_DISK)
        self._error_data.extend(self._explicit_refined_models.filter(lambda seed: seed.error_state < 0).map(self.get_error_metadata).collect())
        self._explicit_refined_models = self._explicit_refined_models.filter(lambda seed: seed.error_state == 0).persist(pyspark.StorageLevel.MEMORY_AND_DISK)

    def retrieve_models(self):
        """
        Retrieve the completed models from the Spark cluster
        :return:
        """
        if not self._constrained_memory:
            finished_models = self.write_models(self._explicit_refined_models.collect())
            return
        n_models = self._explicit_refined_models.count()
        #divide into slices:
        n_full_slices = n_models / self._n_per_slice
        remainder = n_models % self._n_per_slice
        n_per_slice = self._n_per_slice
        for i in range(n_full_slices):
            self.write_models(self._explicit_refined_models.filter(lambda seed: i*n_per_slice <= seed.model_id < (i+1)*n_per_slice).collect())
        self.write_models(self._explicit_refined_models.filter(lambda seed: seed.model_id > n_full_slices*n_per_slice).collect())


    def parallel_write_models(self):
        """
        Call to avoid out of memory errors on driver
        This will cause each worker to write its respective models
        :return:
        """
        models_directory = self._models_directory
        def parallel_map(x):
            SparkDriver.map_write_model(x, models_directory)
            return 0

        collected = self._explicit_refined_models.map(parallel_map).collect()


    def write_models(self, models_to_write):
        """
        Write out the resulting models as system, state, integrator
        :param models_to_write:
        :return:
        """
        os.chdir(self._models_directory)
        for model_seed in models_to_write:
            model_path = os.path.abspath(os.path.join(self._models_directory, model_seed.template_id))
            if not os.path.exists(model_path):
                os.mkdir(model_path)
            try:
                os.chdir(model_path)
                system_file = open('system.xml.gz', mode='w')
                integrator_file = open('integrator.xml.gz', mode='w')
                state_file = open('state.xml.gz', mode='w')
                pdb_file = open('model.pdb', mode='w')
                system_file.writelines(model_seed.explicit_refined_system)
                integrator_file.writelines(model_seed.explicit_refined_integrator)
                state_file.writelines(model_seed.explicit_refined_state)
                pdb_file.writelines(model_seed.explicit_refined_pdb)
            except Exception, e:
                print(str(e))
                continue
            finally:
                pdb_file.close()
                system_file.close()
                integrator_file.close()
                state_file.close()
                os.chdir(self._models_directory)


    def write_model_metadata(self, filename, completed=False):
        """
        Writes out the sequence ID, blast e value, rmsd to reference (if it was computed), and template id

        Parameters
        ----------
        completed : Boolean, default False
             Whether to extract the data from a completed set of models or not
        :return:
        """
        def get_metadata(msmseed):
            result = SparkDriver.get_model_metadata(msmseed)
            return result
        if completed:
            metadata = self._explicit_refined_models.map(get_metadata).collect()
        else:
            metadata = self._modeled_seeds.map(get_metadata).collect()
        header = "TemplateID,blast_eval,seqid,rmsd_to_reference,template_sequence\n"
        metadata_string = "".join(metadata)
        outfile = open(filename, 'w')
        outfile.writelines(metadata_string)
        outfile.close()


    def write_template_models(self):
        """
        maps the function to write out template models for later analysis
        :return:
        """
        template_model_directory = self._template_model_directory
        self._modeled_seeds.map(lambda x: SparkDriver.write_template_model(x, template_model_directory))
        return 0

    @staticmethod
    def get_model_metadata(msmseed):
        try:
            seqid = msmseed.sequence_similarity
        except:
            seqid = "NotComputed"
        blast_eval = str(msmseed.blast_eval)
        template_id = msmseed.template_id
        template_sequence = msmseed.template_sequence
        try:
            rmsd_to_reference = msmseed.rmsd_to_reference
        except:
            rmsd_to_reference = "NotComputed"

        return "%s,%s,%s,%s,%s" % (template_id, blast_eval, seqid, rmsd_to_reference, template_sequence)


    @staticmethod
    def map_write_model(model_seed, model_directory):
        """
        Write out the resulting models as system, state, integrator
        :param models_to_write:
        :return:
        """
        os.chdir(model_directory)
        model_path = os.path.abspath(os.path.join(model_directory, model_seed.template_id))
        if not os.path.exists(model_path):
            os.mkdir(model_path)
        try:
            os.chdir(model_path)
            system_file = open('system.xml.gz', mode='w')
            integrator_file = open('integrator.xml.gz', mode='w')
            state_file = open('state.xml.gz', mode='w')
            pdb_file = open('model.pdb', mode='w')
            system_file.writelines(model_seed.explicit_refined_system)
            integrator_file.writelines(model_seed.explicit_refined_integrator)
            state_file.writelines(model_seed.explicit_refined_state)
            pdb_file.writelines(model_seed.explicit_refined_pdb)
        except Exception, e:
            print(str(e))
        finally:
            pdb_file.close()
            system_file.close()
            integrator_file.close()
            state_file.close()
            os.chdir(model_directory)

    @staticmethod
    def write_template_model(msmseed, template_model_folder):
        """
        This is a utility function for writing out the template models for later examination/reproduction

        Parameters
        ----------
        msmseed : MSMSeed object
            Contains the template model
        template_model_folder : string
            Location of the folder for template models
        """
        import os
        import simtk.openmm.app as app
        os.chdir(template_model_folder)
        template_id = msmseed.template_id
        template_pdb = msmseed.template_structure
        try:
            outfile = open(template_id+'.pdb','w')
            app.PDBFile.writeHeader(template_pdb.topology, file=outfile)
            app.PDBFile.writeModel(template_pdb.topology, template_pdb.positions, file=outfile)
            app.PDBFile.writeFooter(template_pdb.topology, file=outfile)
            outfile.close()
            return msmseed
        except Exception, e:
            msmseed.error_message = str(e)
            msmseed.error_state = -1
            return msmseed


    def write_error_data(self, filename):
        try:
            error_file = open(filename, 'w')
            for error in self._error_data:
                error_file.write(str(error))
            error_file.close()
        except Exception, e:
            print(str(e))

    @staticmethod
    def get_error_metadata(model):
        return [model.template_id, model.error_state, model.error_message]


















