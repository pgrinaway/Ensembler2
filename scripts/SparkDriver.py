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

    def __init__(self, target_squence, spark_master, models_directory, constrained_memory=True, executor_memory='5g', auto_initialize=True,platform = 'CUDA', n_per_slice=10):
        self._target_sequence = target_squence
        spark_conf = pyspark.SparkConf()
        spark_conf.setMaster(spark_master)
        spark_conf.set("spark.executor.memory", executor_memory)
        self._spark_context = pyspark.SparkContext(conf=spark_conf)
        self._models_directory = models_directory
        self._model_rdds = dict()
        self._platform = platform
        self._n_per_slice = n_per_slice
        self._constrained_memory = constrained_memory
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

        self._modeled_seeds = self._rdd_starting_templates.map(modelling.align_template_target).map(modelling.make_model)
        self._failed_seeds = self._modeled_seeds.filter(lambda seed: seed.error_code < 0) #TODO: do something to collect failures
        self._modeled_seeds = self._modeled_seeds.filter(lambda seed: seed.error_code == 0)

    def refine_implicit(self):
        """
        Perform implicit-solvent refinement of the models
        :return:
        """
        self._implicit_refined_seeds = self._modeled_seeds.map(lambda seed: refinement.refine_implicitMD(seed, self._platform, niterations=500, nsteps_per_iteration=100)).persist(pyspark.StorageLevel.MEMORY_AND_DISK)
        self._failed_seeds = self._implicit_refined_seeds.filter(lambda seed: seed.error_code < 0)
        self._implicit_refined_seeds = self._implicit_refined_seeds.filter(lambda seed: seed.error_code == 0).persist(pyspark.StorageLevel.MEMORY_AND_DISK)

    def solvate_models(self):
        """
        Determine the number of waters needed to minimally solvate the model
        :return:
        """
        nwaters_array = self._implicit_refined_seeds.map(solvation.solvate_models).collect()
        nwaters_target = solvation.calculate_nwaters(nwaters_array)
        self._solvated_models = self._implicit_refined_seeds.map(lambda seed: solvation.solvate_models_to_target(seed, nwaters_target))
        self._failed_seeds = self._solvated_models.filter(lambda seed: seed.error_code < 0)
        self._solvated_models = self._solvated_models.filter(lambda seed: seed.error_code == 0).persist(pyspark.StorageLevel.MEMORY_AND_DISK)

    def explicit_refine_models(self):
        """
        Perform explicit refinement of models
        :return:
        """
        self._explicit_refined_models = self._solvated_models.map(lambda seed: refinement.refine_explicitMD(seed, openmm_platform=self._platform, niterations=500, nsteps_per_iteration=100)).persist(pyspark.StorageLevel.MEMORY_AND_DISK)
        self._failed_seeds = self._explicit_refined_models.filter(lambda seed: seed.error_code < 0)
        self._explicit_refined_models = self._explicit_refined_models.filter(lambda seed: seed.error_code == 0).persist(pyspark.StorageLevel.MEMORY_AND_DISK)

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
        for i in range(n_full_slices):
            self.write_models(self._explicit_refined_models.filter(lambda seed: i*self._n_per_slice <= seed.model_id < (i+1)*self._n_per_slice).collect())
        self.write_models(self._explicit_refined_models.filter(lambda seed: seed.model_id > n_full_slices*self._n_per_slice))



    def write_models(self, models_to_write):
        """
        Write out the resulting models as system, state, integrator
        :param models_to_write:
        :return:
        """
        os.chdir(self._models_directory)
        for model_seed in models_to_write:
            model_path = os.path.join(self._models_directory, model_seed.template_id+'_'+model_seed.clone_id)
            if not os.path.exists(model_path):
                os.mkdir(model_path)
            try:
                os.chdir(model_path)
                system_file = open('system.xml.gz', mode='w')
                integrator_file = open('integrator.xml.gz', mode='w')
                state_file = open('state.xml.gz', mode='w')
                pdb_file = open('model.pdb', mode='w')
                system_file.writelines(model_seed.serialized_system)
                integrator_file.writelines(model_seed.serialized_integrator)
                state_file.writelines(model_seed.serialized_state)
                app.PDBFile.writeHeader(model_seed.explicit_refined_pdb.topology, file=pdb_file)
                app.PDBFile.writeModel(model_seed.explicit_refined_pdb.topology,model_seed.explicit_refined_pdb.positions, file=pdb_file)
                app.PDBFile.writeFooter(model_seed.explicit_refined_pdb.topology, file=pdb_file)
            except Exception, e:
                print(str(e))
                continue

    def write_model_metadata(self, model):


















