from pymatgen.core.structure import Molecule

from rxnrlx.jaguar.create_inputs import jaguar_input
from rxnrlx.jaguar.read_files import get_energy_from_file, get_mols_from_irc, verify_success
from rxnrlx.common.utils import sec_to_str

import os, random, subprocess

import time


def ts_relax(ts_guess:Molecule, user_parameters:dict, num_tasks:int) -> Molecule:
    # create new folder for inital TS_relaxation
    os.mkdir("./ts_relaxation")
    os.chdir("./ts_relaxation")

    # Set necessary parameters for code functionality
    user_parameters["ip175"] = 2 # creates XYZ files
    user_parameters["molchg"] = ts_guess.charge 
    user_parameters["multip"] = ts_guess.spin_multiplicity 

    # Create input file
    jaguar_input("ts_opt.in", ts_guess, user_parameters)

    # Submit the job and wait
    start_time = time.time()
    job_id = random.randint(10**8, (10**9)-1)
    process = subprocess.Popen(
        f"$SCHRODINGER/jaguar run -jobname ts_relax_{job_id} -PARALLEL {num_tasks} ts_opt.in -W > ts_opt.out", 
        shell=True
    )
    print("\nRunning 1 Transition State Optimization:")
    print(f"1) $SCHRODINGER/jaguar run ts_opt.in -jobname ts_relax_{job_id} -PARALLEL {num_tasks}")
    process.wait()
    duration = time.time() - start_time

    # Print job result
    if not verify_success("ts_opt.out", "ts_opt"):
        print(f"TS Relaxation failed after: {sec_to_str(duration)}")
        raise Exception("TS Relaxation did not converge")
    else:
        print(f"TS Relaxation finished successfully after: {sec_to_str(duration)}")


    # If the process succeeded, open the optimized TS structure
    opt_ts = Molecule.from_file("ts_opt.xyz") # TODO: double check if this is the correct file to read in
    opt_ts.set_charge_and_spin(charge=ts_guess.charge, spin_multiplicity=ts_guess._spin_multiplicity)
    
    # TODO: Check if the process suceeded or failed
    print("TS Relaxation Succeeded \n")

    # Return to main job folder
    os.chdir("./..")

    return opt_ts


def irc(transition_state:Molecule, user_parameters:dict, num_tasks:int) -> tuple[Molecule, Molecule]:
    # create new folder for inital TS_relaxation
    os.mkdir("./irc_calculation")
    os.chdir("./irc_calculation")

    # Set necessary parameters for code functionality
    user_parameters["babel"] = "xyz" # creates XYZ files
    user_parameters["molchg"] = transition_state.charge 
    user_parameters["multip"] = transition_state.spin_multiplicity 
    
    # Create input file
    jaguar_input("irc.in", transition_state, user_parameters)

    # Submit the job and wait
    start_time = time.time()
    job_id = random.randint(10**8, (10**9)-1)
    process = subprocess.Popen(
        f"$SCHRODINGER/jaguar run -jobname irc_{job_id} -PARALLEL {num_tasks} irc.in -W > irc.out", 
        shell=True
    )
    print("Running 1 IRC Job:")
    print(f"1) $SCHRODINGER/jaguar run irc.in -jobname irc_{job_id} -PARALLEL {num_tasks}")
    process.wait()
    duration = time.time() - start_time

    # Print job result
    if not verify_success("irc.out", "irc"):
        print(f"IRC Calculation failed after: {sec_to_str(duration)}")
        raise Exception("IRC Calculation did not converge")
    else:
        print(f"IRC Calculation finished successfully after: {sec_to_str(duration)}")

    forward_molecule, reverse_molecule = get_mols_from_irc(
        outfile="irc.out", 
        num_atoms=len(transition_state)
    )

    # Add charge and multiplicity information
    forward_molecule.set_charge_and_spin(charge=transition_state.charge, spin_multiplicity=transition_state._spin_multiplicity)
    reverse_molecule.set_charge_and_spin(charge=transition_state.charge, spin_multiplicity=transition_state._spin_multiplicity)

    # Save structures to assist with manual debugging
    forward_molecule.to("forward_molecule.xyz")
    reverse_molecule.to("reverse_molecule.xyz")

    # Return to main job folder
    os.chdir("./..")

    return forward_molecule, reverse_molecule

def geom_opt(forward_molecule, reverse_molecule, user_parameters, num_tasks) -> tuple[Molecule, Molecule]:

    # create new folder for inital TS_relaxation
    os.mkdir("./geometry_optimizations")
    os.chdir("./geometry_optimizations")

    # Set necessary parameters for code functionality
    user_parameters["ip175"] = 2 # creates XYZ files

    # Run a geometry opt for each molecule
    print("Running 2 Optimizations:")
    processes = list()
    start_time = time.time()
    for i, (molec, ext) in enumerate(zip([forward_molecule, reverse_molecule], ["fwd", "rev"])):
        # Set charge and multiplicity
        user_parameters["molchg"] = molec.charge
        user_parameters["multip"] = molec.spin_multiplicity

        # Create input file
        jaguar_input(f"opt_{ext}.in", molec, user_parameters)

        job_id = random.randint(10**8, (10**9)-1)
        subproc = subprocess.Popen(
            f"$SCHRODINGER/jaguar run -jobname opt_{ext}_{job_id} -PARALLEL {num_tasks//2} opt_{ext}.in -W > opt_{ext}.out", 
            shell=True
        )
        print(f"{i+1}) $SCHRODINGER/jaguar run opt_{ext}.in -jobname opt_{ext}_{job_id} -PARALLEL {num_tasks//2}")
        processes.append(subproc)
    
    for subp in processes:
        subp.wait()

    duration = time.time() - start_time

    # Print overall job result
    print(f"Optimization jobs finished after: {sec_to_str(duration)}")

    # Print more specific job results
    fwd_result = verify_success("opt_fwd.out", "opt_fwd")
    rev_result = verify_success("opt_rev.out", "opt_rev")
    
    print(f"Forward Molecule Optimiation: {'SUCCESSFUL' if fwd_result else 'FAILED'}")
    print(f"Reverse Molecule Optimization: {'SUCCESSFUL' if rev_result else 'FAILED'}")
    
    if not (fwd_result and rev_result):
        raise Exception("At least one geometry optimization failed")



    # Read the structures into Molecule objects
    fwd = Molecule.from_file("opt_fwd.xyz") # TODO: double check if this is the correct file to read in
    fwd.set_charge_and_spin(charge=forward_molecule.charge, spin_multiplicity=forward_molecule._spin_multiplicity)
    rev = Molecule.from_file("opt_rev.xyz") # TODO: double check if this is the correct file to read in
    fwd.set_charge_and_spin(charge=reverse_molecule.charge, spin_multiplicity=reverse_molecule._spin_multiplicity)

    print("Geometry Optimizations Finished\n")

    # Return to folder
    os.chdir("./..")

    return fwd, rev


def calculate_gibbs(forward_molecule, reverse_molecule, transition_state, user_parameters, num_tasks) -> dict:
    """ Submits a single point energy calculation """

    os.mkdir("./energy_calculation")
    os.chdir("./energy_calculation")

    # Run a single point calculation for each molecule
    print("Running 3 Frequency Calculations")
    processes = list()
    for i, (molec, ext) in enumerate(zip([forward_molecule, reverse_molecule, transition_state], ["fwd", "rev", "ts"])):
        # set charge and multiplicity
        user_parameters["molchg"] = molec.charge
        user_parameters["multip"] = molec.spin_multiplicity

        # Create input file
        jaguar_input(f"energy_{ext}.in", molec, user_parameters)

        job_id = random.randint(10**8, (10**9)-1)
        subproc = subprocess.Popen(
            f"$SCHRODINGER/jaguar run -jobname energy_{ext}_{job_id} -PARALLEL {num_tasks//3} energy_{ext}.in -W > energy_{ext}.out", 
            shell=True
        )
        print(f"{i+1}) $SCHRODINGER/jaguar run energy_{ext}.in -jobname energy_{ext}_{job_id} -PARALLEL {num_tasks//3}")
        processes.append(subproc)
    
    for subp in processes:
        subp.wait()
    
    fwd_result = verify_success("energy_fwd.out", "energy_fwd")
    rev_result = verify_success("energy_rev.out", "energy_rev")
    ts_result = verify_success("energy_ts.out", "energy_ts")

    print(f"Forward Molecule Freqency Calculation: {'SUCCESSFUL' if fwd_result else 'FAILED'}")
    print(f"Reverse Molecule Freqency Calculation: {'SUCCESSFUL' if rev_result else 'FAILED'}")
    print(f"Transition State Freqency Calculation: {'SUCCESSFUL' if ts_result else 'FAILED'}")

    if not (fwd_result and rev_result and ts_result):
        raise Exception("At least one frequency calculation failed")

    # Get energetics from each outfile
    forward_energy = get_energy_from_file("energy_fwd.out")
    reverse_energy = get_energy_from_file("energy_rev.out")
    ts_energy = get_energy_from_file("energy_ts.out")

    # Return to main folder
    os.chdir("..")

    return {
        "forward": forward_energy,
        "reverse": reverse_energy,
        "transition_state": ts_energy
    }