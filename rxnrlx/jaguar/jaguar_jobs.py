from pymatgen.core.structure import Molecule

from rxnrlx.jaguar.create_inputs import jaguar_input

import os, random, subprocess, re


def ts_relax(ts_guess:Molecule, user_parameters:dict, num_tasks:int) -> Molecule:
    # create new folder for inital TS_relaxation
    os.mkdir("./initial_ts_relaxation")
    os.chdir("./initial_ts_relaxation")
    
    # create TS relaxation input file
    ts_parameters = {
        "igeopt" : 2,
        "inhess": 4,
        "epsout": 18.5,
        "isolv": 7,
        "maxitg": 300,
        "isymm": 0,
        "basis": "DEF2-SVPD",
        "maxit": 300,
        "ip472": 2,
        "dftname": "wB97X-D",
        "nogas": 2,
        "ip175": 2,
        "iacc": 2,
        "molchg": ts_guess.charge,
        "multip": ts_guess.spin_multiplicity
    }
    for key, val in user_parameters: # Update with any user-specified parameters
        ts_parameters[key] = val
    jaguar_input("ts_opt.in", ts_guess, ts_parameters)

    # Submit the job and wait
    job_id = random.randint(10^8, 10^9-1)
    process = subprocess.Popen(["$SCHRODINGER/jaguar", f"run -jobname ts_relax_{job_id} -PARALLEL {num_tasks} ts_opt.in -W"])
    print("Initial TS Relaxation Started")
    process.wait()

    # TODO: Check if the process suceeded or failed
    print("TS Relaxation Finished\n")

    # If the process succeeded, open the optimized TS structure
    opt_ts = Molecule.from_file("ts_opt.xyz") # TODO: double check if this is the correct file to read in
    opt_ts.set_charge_and_spin(charge=ts_guess.charge, spin_multiplicity=ts_guess._spin_multiplicity)
    
    # Return to main job folder
    os.chdir("/..")

    return opt_ts


def irc(transition_state:Molecule, user_parameters:dict, num_tasks:int) -> tuple[Molecule, Molecule]:
    # create new folder for inital TS_relaxation
    os.mkdir("./irc_calculation")
    os.chdir("./irc_calculation")

    irc_parameters = {
        "inhess": 4,
        "no_mul_imag_freq": 1,
        "epsout": 18.5,
        "valid_sections": 0,
        "isolv": 7,
        "maxitg": 60000,
        "isymm": 0,
        "basis": "def2-svpd",
        "ircstep": 0.1,
        "maxit": 300,
        "geoconv_mode": "standard",
        "itrvec": 1,
        "babel": "xyz",
        "ircmxcyc": 300,
        "dftname": "wb97x-d",
        "lqa_step": 1,
        "nogas": 2,
        "scale_geoconv": 3.0,
        "ircmax": 100,
        "irc": 1,
        "molchg": transition_state.charge,
        "multip": transition_state.spin_multiplicity
    }
    for key, val in user_parameters: # Update with any user-specified parameters
        irc_parameters[key] = val
    jaguar_input("irc.in", transition_state, irc_parameters)

    # Submit the job and wait
    job_id = random.randint(10^8, 10^9-1)
    process = subprocess.Popen(["$SCHRODINGER/jaguar", f"run -jobname irc_{job_id} -PARALLEL {num_tasks} irc.in -W"])
    print("IRC Calculation Started")
    process.wait()
    print("IRC Calculation Finished\n")


    forward_molecule, reverse_molecule = get_mols_from_irc(
        outfile="irc.out", 
        num_atoms=len(transition_state)
    )

    # Add charge and multiplicity information
    forward_molecule.set_charge_and_spin(charge=transition_state.charge, spin_multiplicity=transition_state._spin_multiplicity)
    reverse_molecule.set_charge_and_spin(charge=transition_state.charge, spin_multiplicity=transition_state._spin_multiplicity)

    # Return to main job folder
    os.chdir("/..")

    return forward_molecule, reverse_molecule

def geom_opt(forward_molecule, reverse_molecule, user_parameters, num_tasks):

    # create new folder for inital TS_relaxation
    os.mkdir("./geometry_optimizations")
    os.chdir("./geometry_optimizations")

    opt_parameters = {
        "igeopt": 1,
        "isolv": 7,
        "isymm": 0,
        "epsout": 18.5,
        "dftname": "wb97x-v",
        "basis": "def2-svpd",
        "babel": "xyz",
        "maxit": 300,
        "maxitg": 300,
        "iacc": 2,
        "nogas": 2,
        "ip175": 2,
        "ip142": 2,
    }
    for key, val in user_parameters: # Update with any user-specified parameters
        opt_parameters[key] = val


    # Run a geometry opt for each molecule
    processes = list()
    for molec, ext in zip([forward_molecule, reverse_molecule], ["fwd", "rev"]):
        opt_parameters["molchg"] = molec.charge
        opt_parameters["multip"] = molec.spin_multiplicity

        jaguar_input(f"opt_{ext}.in", molec, opt_parameters)

        job_id = random.randint(10^8, 10^9-1)
        subproc = subprocess.Popen(["$SCHRODINGER/jaguar", f"run -jobname irc_{job_id} -PARALLEL {num_tasks//2} opt_{ext}.in -W"])
        processes.append(subproc)
    
    # Wait for both to finish
    print("Geometry Optimizations Started")
    for subp in processes:
        subp.wait()
    print("Geometry Optimizations Finished\n")

    # Read the structures into Molecule objects
    fwd = Molecule.from_file("opt_fwd.xyz") # TODO: double check if this is the correct file to read in
    fwd.set_charge_and_spin(charge=forward_molecule.charge, spin_multiplicity=forward_molecule._spin_multiplicity)
    rev = Molecule.from_file("opt_rev.xyz") # TODO: double check if this is the correct file to read in
    fwd.set_charge_and_spin(charge=reverse_molecule.charge, spin_multiplicity=reverse_molecule._spin_multiplicity)

    # Return to folder
    os.chdir("/..")

    return fwd, rev

def get_mols_from_irc(outfile, num_atoms) -> tuple[Molecule, Molecule]:
    """ Get the optimized forward and backward molecules from the transition state """
    
    with open(outfile, "r") as f:
        lines = f.readlines()

    # Get the places to search for the geometry definitions
    for i, line in enumerate(lines):
        if "Forward IRC cycle complete" in line:
            forward_section = i
        if "Reverse IRC cycle complete" in line:
            reverse_section = i
            break

    # find the molecules
    forward_molecule = find_molecule_in_section(lines, forward_section, num_atoms)
    reverse_molecule = find_molecule_in_section(lines, reverse_section, num_atoms)

    return forward_molecule, reverse_molecule


def find_molecule_in_section(lines, starting_place, num_atoms) -> Molecule:
    """ Find the first relaxed molecule definition to appear after the given line index """

    # Find the geometry header line
    found = False
    i = starting_place
    while found:
        pattern = re.compile(r"atom\s+x\s+y\s+z")
        if re.match(pattern, lines[i]):
            found = True
        i += 1

   
    # i should now hold the first line of the atoms
    species_list = list()
    coord_list = list()
    for j in range(num_atoms):
        species, x, y, z = lines[i+j].split()

        # Remove species index
        species = re.sub(r'[^a-zA-Z]', '', species)

        # add values to lists
        species_list.append(species)
        coord_list.append([float(x), float(y), float(z)])

    # Return read in molecule (no charge or multiplicity information)
    return Molecule(
        species=species_list,
        coords=coord_list,
    )