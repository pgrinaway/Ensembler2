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
    import simtk.openmm as openmm
    import simtk.unit as unit
    import gzip
    import StringIO

    #get the serialized system, state, integrator
    serialized_system = "".join(gzip.GzipFile(fileobj=StringIO.StringIO(msmseed.explicit_refined_system)).readlines())
    serialized_state = "".join(gzip.GzipFile(fileobj=StringIO.StringIO(msmseed.explicit_refined_state)).readlines())
    serialized_integrator = "".join(gzip.GzipFile(fileobj=StringIO.StringIO(msmseed.explicit_refined_integrator)).readlines())
    system = openmm.XmlSerializer.deserialize(serialized_system)
    #integrator = openmm.XmlSerializer.deserialize(serialized_integrator)
    state = openmm.XmlSerializer.deserialize(serialized_state)

    #integrator parameters
    timestep = 2.0 * unit.femtoseconds # timestep
    temperature = 300.0 * unit.kelvin # simulation temperature
    pressure = 1.0 * unit.atmospheres # simulation pressure
    collision_rate = 20.0 / unit.picoseconds # Langevin collision rate
    barostat_period = 50

    #set the temperatures + make new MSMSeed objects
    clone_list = list()
    for i in range(n_clones):
        integrator = openmm.LangevinIntegrator(temperature, collision_rate, timestep)
        context = openmm.Context(system, integrator)
        context.setState(state)
        context.setVelocitiesToTemperature(temperature)
        new_state = context.getState(getPositions=True, getVelocities=True, getForces=True, getEnergy=True, getParameters=True, enforcePeriodicBox=True)


    #write out new serialized system, state, integrator