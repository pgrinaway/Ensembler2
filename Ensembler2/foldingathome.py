__author__ = 'Patrick B. Grinaway'

def package_for_fah(msmseed, n_clones, retain_metadata = False):
    """
    This generates clones for a folding at home run, each of which have different randomized velocities

    Arguments
    ---------
    msmseed : an MSMSeed object
        Contains the explicit-solvent refined model to use
    n_clones : int
        the number of clones to generate for each model
    retain_metadata : Boolean, optional
        whether to transfer all intermediate and metadata to the clone objects.
        default false to save memory

    Returns
    -------
    fah_clone : a n_clones long list of MSMSeed objects
         MSMSeed objects with serialized state, integrator, system, and pdb model

    """
    raise NotImplementedError