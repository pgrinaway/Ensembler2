__author__ = 'Patrick B. Grinaway'



if __name__=="__main__":
    import Ensembler2.SparkDriver as SparkDriver
    f = open('target_seq.fasta', 'r')
    tgt_seq = "".join(f.readlines())
    driver = SparkDriver.SparkDriver(tgt_seq, 'spark://gpu-2-13.local:7077', 'output_models', executor_memory='5g', auto_initialize=False)
    driver.initialize_with_blast(max_hits=10)
    driver.make_models()
    driver.refine_implicit()
    driver.solvate_models()
    driver.explicit_refine_models()
    driver.retrieve_models()
    driver.write_error_data('errored_models.txt')