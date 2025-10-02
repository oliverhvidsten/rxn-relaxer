"""
TODO:
The options for this file should be as follow:
- Reoptimize the molecules with a new functional/basis set
- Calculate Gibbs Free Energies
"""
from pymatgen.core.structure import Molecule

from rxnrlx.common.constants import FWD_FILENAME, REV_FILENAME, TS_FILENAME

import os, sys, yaml

def refine(config):
    """
    The options for this file should be as follow:
    - Reoptimize the molecules with a new functional/basis set
    - Calculate Gibbs Free Energies
    """ 

    # User has the option to specify the old job folder or individual molecules
    # If they specify the old_job_folder, this program will grab the species from the final_structures subfolder
    if "old_job_folder" in config["info"]:
        os.chdir(config["info"]["old_job_folder"])
        forward_molecule = Molecule.from_file(f'./final_structures/{FWD_FILENAME}')
        reverse_molecule = Molecule.from_file(f'./final_structures/{REV_FILENAME}')
        transition_state = Molecule.from_file(f'./final_structures/{TS_FILENAME}')

    else: # In the event they did not specify the job folder, they should have specified 3 file locations where the molecules are located
        if ("forward" in config["info"]) and ("reverse" in config["info"]) and ("transition_state" in config["info"]):
            forward_molecule = Molecule.from_file(config["info"]["forward"])
            reverse_molecule = Molecule.from_file(config["info"]["reverse"])
            transition_state = Molecule.from_file(config["info"]["transition_state"])

        else:
            raise Exception("The info section of the config file should contain either {\'old_job_folder\'} or {\'forward\', \'reverse\', and \'transition_state\'}")
    
    # Set charge and multiplicity values for the molecules
    for molec in [forward_molecule, reverse_molecule, transition_state]:
        molec.set_charge_and_spin(
            charge=config["info"].get("charge", 0),
            spin_multiplicity=config["info"].get("multiplicity", 1)
        )

    # implementation
    if config["info"]["software"] == "jaguar": 
        from rxnrlx.jaguar.jaguar_jobs import ts_relax, geom_opt, calculate_gibbs
    else:
        raise NotImplementedError()
    
    os.mkdir("./refine_structures")
    os.chdir("./refine_structures")

    # If user requests re-optimization of the inputs:
    if config["info"]["reoptimize"]:
        stable_refined = False
        ts_refined = False
        try:
            # Relax TS with new level of theory
            transition_state = ts_relax(
                ts_guess=transition_state, 
                user_parameters=config.get("ts_relax", {}),
                num_tasks=config["info"].get("ntasks", 2)
            )
        except:
            if config["info"].get("die_on_ts_failure", True):
                raise Exception("TS Optimization Failed")
            else:
                print("TS Optimization failed, keeping original geometry for energy calculation")
        else:
            ts_refined = True # New TS was optimized with this level of theory
        
        # Optimize stable reactant and product with new level of theory 
        try:
            forward_molecule, reverse_molecule = geom_opt(
                forward_molecule=forward_molecule,
                reverse_molecule=reverse_molecule,
                user_parameters=config.get("geom_opt", {}),
                num_tasks=config["info"].get("ntasks", 2)
            )
        except:
            if config["info"].get("die_on_ts_failure", True):
                raise Exception("Product and/or Reactant Optimizations Failed")
            else:
                print("Stable Geometry Optimizations Failed, keeping original geometries for energy calculation")
        else:
            stable_refined = True # New stable geometries optimized with this level of theory

        ## Save refined structures
        os.mkdir("./refined_final_structures")
        os.chdir("./refined_final_structures")
        if stable_refined:
            forward_molecule.to(FWD_FILENAME)
            reverse_molecule.to(REV_FILENAME)
        if ts_refined:
            transition_state.to(TS_FILENAME)

    # Next run the energetic calculations
    energy_info = calculate_gibbs(
        forward_molecule=forward_molecule, 
        reverse_molecule=reverse_molecule, 
        transition_state=transition_state, 
        user_parameters=config.get("energy"), 
        num_tasks=config["info"].get("ntasks")
        )
    
    # Get reaction energetic information in electron Volts (eV)
    dG = (energy_info["forward"] - energy_info["reverse"]) * 27.114 
    barrier = (energy_info["transition_state"] - energy_info["reverse"]) * 27.114 
    reverse_barrier = (energy_info["transition_state"] - energy_info["forward"]) * 27.114 

    energy_info["Reaction Info (eV)"] = {
        "Delta G": dG, 
        "Forward Activation Barrier": barrier,
        "Reverse Activation Barrier": reverse_barrier
        }



    with open("energy.yaml", 'w') as f:
        yaml.dump(energy_info, f, default_flow_style=False)


if __name__ == "__main__":
    """ Read in Command Line Arguments """
    if len(sys.argv) == 2:
        print(f"Configuration File: {sys.argv[1]}")
        try:
            with open(sys.argv[1], "r") as f:
                config = yaml.safe_load(f)
        except:
            raise Exception("Invalid File Specified")
    else:
        error_message = [f"Invalid number of arguments ({len(sys.argv)}).",
                         "Use format: ts2rxn.py <config_file.yaml>"]
        raise Exception("\n".join(error_message))
        
    # Run main code
    refine(config)