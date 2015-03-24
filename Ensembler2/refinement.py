__author__ = 'Patrick B. Grinaway'


def refine_implicitMD(msmseed, openmm_platform='CPU', niterations=100, nsteps_per_iteration=5):
    """
    Use OpenMM to perform an implicit solvent calculation on the output of MODELLER.

    Parameters
    ----------
    msmseed : MSMSeed
        object containing the modeled structure to be used for implicit simulation
    openmm_platofm : String, optional
        name of the OpenMM platform to be used in the simulation. Default CPU.
    niterations : int, optional
        number of iterations of steps to perform when running dynamics. Default 100.
    nsteps_per_iteration : int, optional
        number of steps to take each iteration. Default 5



    """
    import os
    import simtk.openmm as openmm
    import simtk.unit as unit
    import simtk.openmm.app as app
    from Ensembler2.MSMSeed import MDSys

    forcefields_to_use = ['amber99sbildn.xml', 'amber99_obc.xml'] # list of forcefields to use in parameterization

    timestep = 2.0 * unit.femtoseconds # timestep
    temperature = 300.0 * unit.kelvin # simulation temperature
    collision_rate = 20.0 / unit.picoseconds # Langevin collision rate
    cutoff = None # nonbonded cutoff
    minimization_tolerance = 10.0 * unit.kilojoules_per_mole / unit.nanometer
    minimization_steps = 20

    kB = unit.MOLAR_GAS_CONSTANT_R
    kT = kB * temperature

    pH = 8.0

    forcefield = app.ForceField(*forcefields_to_use)
    platform = openmm.Platform.getPlatformByName(openmm_platform)
        #this will just be None if there is no gpu
    gpuid = os.getenv("CUDA_VISIBLE_DEVICES")
    if openmm_platform == 'CUDA':
        #here, gpuid is returned as a string, so it can be directly given to the platform method
        platform.setPropertyDefaultValue('CudaDeviceIndex', gpuid)
    if openmm_platform == 'OpenCL':
        platform.setPropertyDefaultValue('OpenCLDeviceIndex', gpuid)
    try:
        modeller = app.Modeller(msmseed.target_model.topology, msmseed.target_model.positions)
        modeller.addHydrogens(forcefield, pH=pH, variants=None)
        topology = modeller.getTopology()
        positions = modeller.getPositions()

        system = forcefield.createSystem(topology, nonbondedMethod=app.NoCutoff, constraints=app.HBonds)
    except Exception, e:
        msmseed.error_state = -2
        msmseed.error_message = str(e)
        return msmseed

    integrator = openmm.LangevinIntegrator(temperature, collision_rate, timestep)
    context = openmm.Context(system, integrator, platform)
    context.setPositions(positions)

    openmm.LocalEnergyMinimizer.minimize(context, minimization_tolerance, minimization_steps)

    for iteration in range(niterations):
            # integrate dynamics
            integrator.step(nsteps_per_iteration)
            # get current state
            state = context.getState(getEnergy=True, getPositions=True)
            simulation_time = state.getTime()
            potential_energy = state.getPotentialEnergy()
            kinetic_energy = state.getKineticEnergy()
            import numpy
            if numpy.isnan(potential_energy/kT) or numpy.isnan(kinetic_energy/kT):
                msmseed.error_message = "Potential or kinetic energies are nan."
                msmseed.error_state=-2
                return msmseed
    state = context.getState(getPositions=True)
    refined = MDSys(modeller.topology, state.getPositions())
    msmseed.implicit_refined_model = refined
    return msmseed

def refine_explicitMD(msmseed, openmm_platform='CPU', niterations=1, nsteps_per_iteration=5):
    """
    Run an explicit-solvent MD refinement on a solvated model structure.

    Parameters
    ----------
    msmseed : MSMSeed
        object containing solvated model to simulate
    openmm_platform : String, optional
        name of the OpenMM platform to use in the simulation. Default CPU
    niterations : int, optional
        number of iterations of integrator steps to take. Default 1.
    nsteps_per_iteration : int, optional
        number of steps to take per iteration. Default 5.

    Returns
    -------
    msmseed : MSMSeed
         object containing gzipped explicitly-refined model pdb, integrator xml, system xml, and state xml

    """
    import simtk.openmm as openmm
    import simtk.openmm.app as app
    import simtk.unit as unit
    import time
    import StringIO
    import gzip
    import os
    import numpy
    from Ensembler2.MSMSeed import MDSys

    platform = openmm.Platform.getPlatformByName(openmm_platform)
    #this will just be None if there is no gpu
    gpuid = os.getenv("CUDA_VISIBLE_DEVICES")
    if openmm_platform == 'CUDA':
        #here, gpuid is returned as a string, so it can be directly given to the platform method
        platform.setPropertyDefaultValue('CudaDeviceIndex', gpuid)
    if openmm_platform == 'OpenCL':
        platform.setPropertyDefaultValue('OpenCLDeviceIndex', gpuid)

    forcefields_to_use = ['amber99sbildn.xml', 'tip3p.xml'] # list of forcefields to use in parameterization

    timestep = 2.0 * unit.femtoseconds # timestep
    temperature = 300.0 * unit.kelvin # simulation temperature
    pressure = 1.0 * unit.atmospheres # simulation pressure
    collision_rate = 20.0 / unit.picoseconds # Langevin collision rate
    barostat_period = 50
    niterations = 100 # number of iterations

    nonbondedMethod = app.PME

    minimization_tolerance = 10.0 * unit.kilojoules_per_mole / unit.nanometer
    minimization_steps = 20

    kB = unit.MOLAR_GAS_CONSTANT_R
    kT = kB * temperature


    forcefield = app.ForceField(*forcefields_to_use)
    solvated_model = msmseed.solvated_model
    system = forcefield.createSystem(solvated_model.topology, nonbondedMethod=nonbondedMethod, constraints=app.HBonds)
    barostat = openmm.MonteCarloBarostat(pressure, temperature, barostat_period)
    system.addForce(barostat)
    integrator = openmm.LangevinIntegrator(temperature, collision_rate, timestep)
    context = openmm.Context(system, integrator, platform)
    context.setPositions(solvated_model.positions)
    openmm.LocalEnergyMinimizer.minimize(context, minimization_tolerance, minimization_steps)


    context.setVelocitiesToTemperature(temperature)
    initial_time = time.time()
    ns_per_day = 0.0
    for iteration in range(niterations):
            # integrate dynamics
            integrator.step(nsteps_per_iteration)
            # get current state
            state = context.getState(getEnergy=True, getPositions=True)
            potential_energy = state.getPotentialEnergy()
            kinetic_energy = state.getKineticEnergy()
            if numpy.isnan(potential_energy/kT) or numpy.isnan(kinetic_energy/kT):
                msmseed.error_status = -5
                msmseed.error_message = 'The simulation has encountered NaNs in energies'
                break
            simulation_time = state.getTime()
            potential_energy = state.getPotentialEnergy()
            kinetic_energy = state.getKineticEnergy()
            final_time = time.time()
            elapsed_time = (final_time - initial_time) * unit.seconds
            ns_per_day = (simulation_time / elapsed_time) / (unit.nanoseconds / unit.day)
            box_vectors = state.getPeriodicBoxVectors()
            volume_in_nm3 = (box_vectors[0][0] * box_vectors[1][1] * box_vectors[2][2]) / (unit.nanometers**3) # TODO: Use full determinant
            remaining_time = elapsed_time * (niterations-iteration-1) / (iteration+1)
    msmseed.explicit_ns_per_day = ns_per_day
    state = context.getState(getPositions=True)
    try:
        #save the pdb of the current model
        explicit_model_pdb = StringIO.StringIO()
        with gzip.GzipFile(fileobj = explicit_model_pdb, mode = 'w') as output:
            app.PDBFile.writeHeader(solvated_model.topology, file=output)
            app.PDBFile.writeModel(solvated_model.topology, state.getPositions(), file=output)
            app.PDBFile.writeFooter(solvated_model.topology, file=output)
        msmseed.explicit_refined_pdb = explicit_model_pdb.getvalue()

    #save the state
        serialized_state_stringio = StringIO.StringIO()
        state_xml = openmm.XmlSerializer.serialize(state)
        with gzip.GzipFile(fileobj=serialized_state_stringio, mode = 'w') as output:
             output.write(state_xml)
        msmseed.explicit_refined_state = serialized_state_stringio.getvalue()

    #save the integrator
        serialized_integrator_stringio = StringIO.StringIO()
        integrator_xml = openmm.XmlSerializer.serialize(integrator)
        with gzip.GzipFile(fileobj = serialized_integrator_stringio, mode = 'w') as buffer:
            buffer.write(integrator_xml)
        msmseed.explicit_refined_integrator = serialized_integrator_stringio.getvalue()

    #save the system
        serialized_system_stringio = StringIO.StringIO()
        system_xml = openmm.XmlSerializer.serialize(system)
        with gzip.GzipFile(fileobj = serialized_system_stringio, mode = 'w') as buffer:
           buffer.write(system_xml)
        msmseed.explicit_refined_system = serialized_system_stringio.getvalue()
    except:
        msmseed.error_state = -5
        msmseed.error_messsage = "Failed to save model"
    finally:
        return msmseed