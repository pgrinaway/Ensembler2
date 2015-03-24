__author__ = 'Patrick B. Grinaway'




def solvate_models(msmseed):
    """
    Calculate the number of waters needed to solvate the implicitly-refined model

    Parameters
    ----------
    msmseed : MSMSeed
        object containing the unsolvated, implicitly refined model structure

    Returns
    -------
    msmseed : MSMSeed
        object additionally containing the number of waters required to solvate the model

    """
    import os
    import simtk.unit as unit
    import simtk.openmm.app as app


    forcefields_to_use = ['amber99sbildn.xml', 'tip3p.xml'] # list of forcefields to use in parameterization
    nparticles_per_water = 3
    padding = 10.0 * unit.angstroms

    forcefield = app.ForceField(*forcefields_to_use)
    topology = msmseed.implicit_refined_model.topology
    positions = msmseed.implicit_refined_model.positions
    natoms_initial = len(positions)
    modeller = app.Modeller(topology, positions)
    modeller.addSolvent(forcefield, model='tip3p', padding=padding)
    positions = modeller.getPositions()
    natoms_final = len(positions)
    nwaters = (natoms_final - natoms_initial) / nparticles_per_water
    msmseed.nwaters = nwaters
    return msmseed



def calculate_nwaters(nwaters_list):
    """
    Calculate the number of waters at the 68th percentile.

    Parameters
    ----------
    nwaters_list : list of ints
        list of the number of waters needed to solvate each model

    Returns
    -------
    index68 : int
        the target number of waters at the 68th percentile
    """
    import numpy as np
    nwaters_array = np.array(nwaters_list)
    nwaters_array.sort()
    index68 = int((len(nwaters_array) - 1) * 0.68)
    #index95 = int((len(nwaters_array) - 1) * 0.95)
    return nwaters_array[index68]

def solvate_models_to_target(msmseed, target_nwaters):
    """
    Solvate the model to the target number of waters. If the model requires more than target_nwaters for solvation,
    it will be rejected.

    Parameters
    ----------
    msmseed : MSMSeed
        object containing an implicitly refined model to solvate
    target_nwaters : int
        number of waters to add

    Returns
    -------
    msmseed : MSMSeed
        object containing a solvated model with target_nwaters, or containing an error message

    """

    import os
    from MSMSeed import MDSys
    import simtk.openmm.app as app
    import simtk.unit as unit
    natoms_per_solvent = 3
    refined_model = msmseed.implicit_refined_model
    # Count initial atoms.
    natoms_initial = len(refined_model.positions)
    forcefields_to_use = ['amber99sbildn.xml', 'tip3p.xml']
    model='tip3p'
    forcefield = app.ForceField(*forcefields_to_use)

    # Solvate with zero padding to determine min number of waters and minimal unit cell dimensions.
    modeller = app.Modeller(refined_model.topology, refined_model.positions)
    modeller.addSolvent(forcefield, model=model, padding=0.0*unit.angstroms)
    topology = modeller.getTopology()
    positions = modeller.getPositions()
    box_min = topology.getUnitCellDimensions()
    natoms_min = len(positions) # minimal number of atoms
    nwaters_min = (natoms_min - natoms_initial) / natoms_per_solvent # minimal number of waters
    volume_min = box_min[0] * box_min[1] * box_min[2]
    residues = [ r for r in topology.residues() ] # build a list of residues
    nresidues_min = len(residues) # number of residues
    if nwaters_min > target_nwaters:
        msmseed.error_state = -3
        msmseed.error_message = "The minimally solvated model has more than the target number of waters"
        return msmseed

    #Make a slightly enlarged box
    scale = 1.1
    modeller = app.Modeller(refined_model.topology, refined_model.positions)
    topology = modeller.getTopology()
    topology.setUnitCellDimensions(box_min * scale)
    modeller.addSolvent(forcefield, model=model)
    positions = modeller.getPositions()
    box_enlarged = topology.getUnitCellDimensions()
    natoms_enlarged = len(positions) # minimal number of atoms
    nwaters_enlarged = (natoms_enlarged - natoms_initial) / natoms_per_solvent # minimal number of waters
    volume_enlarged = box_enlarged[0] * box_enlarged[1] * box_enlarged[2]
    density = (nwaters_enlarged - nwaters_min) / (volume_enlarged - volume_min)
    over_target = False
    extra_nwaters = 100
    while not over_target:
        delta_volume = (target_nwaters + extra_nwaters - nwaters_min) / density
        scale = ((volume_min + delta_volume) / volume_min)**(1.0/3.0)
        modeller = app.Modeller(refined_model.topology, refined_model.positions)
        topology = modeller.getTopology()
        topology.setUnitCellDimensions(box_min * scale)
        modeller.addSolvent(forcefield, model=model)
        positions = modeller.getPositions()
        topology = modeller.getTopology()
        natoms = len(positions) # minimal number of atoms
        nwaters = (natoms - natoms_initial) / natoms_per_solvent # minimal number of waters
        if (nwaters > target_nwaters):
            over_target = True
        else:
            extra_nwaters += 100

    # Delete waters to achieve target.
    ndelete = nwaters - target_nwaters
    if (ndelete > 0):
        residues = [ r for r in topology.residues() ] # build a list of residues
        nresidues = len(residues)

        # Select a random subset to delete.
        import numpy.random
        indices = numpy.random.permutation(range(nresidues_min,nresidues))
        residues_to_delete = list()
        for index in indices[0:ndelete]:
            residues_to_delete.append(residues[index])

        modeller.delete(residues_to_delete)
     # Get topology and positions.
        topology = modeller.getTopology()
        positions = modeller.getPositions()

        # Count number of waters.
        natoms_final = len(positions)
        nwaters = (natoms_final - natoms_initial) / 3
    solvated_mdsys = MDSys(topology, positions)
    msmseed.solvated_model = solvated_mdsys
    return msmseed


