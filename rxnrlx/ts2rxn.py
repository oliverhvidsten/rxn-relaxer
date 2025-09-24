from pymatgen.core.structure import Molecule
import os, sys, yaml

def ts2rxn(config:dict={}):
    """
    This function orchestrates a workflow that takes a ts_guess and turns it into a reaction pathway.
    Steps
    - Optimize TS Guess
    - Perform IRC Analysis on optimized TS
    - Perform geometry optimizations on the forward and reverse IRC structures
    - Output a diagram that shows the energetic pathway of the reaction (reactant, TS, product)
        - assume the reactant is the "reverse" direction from IRC
    """
    # Get filename from config
    ts_guess_file = config["info"]["ts_guess_filename"]
    
    # open xyz file from args
    ts_guess = Molecule.from_file(ts_guess_file)
    ts_guess.set_charge_and_spin(
        charge=config["info"].get("charge", 0),
        spin_multiplicity=config["info"].get("multiplicity", 1)
        )

    # create folder for job to be run in and move into the directory
    job_name = config["info"]["job_name"]
    os.makedirs(f"./{job_name}")
    os.chdir(f"./{job_name}")

    # implementation
    if config["info"]["software"] == "jaguar": 
        from rxnrlx.jaguar.jaguar_jobs import ts_relax, irc, geom_opt
    else:
        raise NotImplementedError()


    # Perform Transition State Optimization
    try:
        transition_state = ts_relax(
            ts_guess=ts_guess, 
            user_parameters=config.get("ts_relax", {}), 
            num_tasks=config.get("submit", {"ntasks":2})["ntasks"]
        )
    except Exception as e:
        # If TS optimization fails, still keep the program going with the guess as the transition state
        if not config["info"].get("die_on_ts_failure", True):
            transition_state = ts_guess
        else:
            raise e


    # Perform IRC Analysis
    forward_molecule, reverse_molecule = irc(
        transition_state=transition_state, 
        user_parameters=config.get("irc", {}), 
        num_tasks=config.get("submit", {"ntasks":2})["ntasks"]
    )


    # Optimize Forward and Reverse Molecule (reactants and products)
    forward_optimized, reverse_optimized = geom_opt(
        forward_molecule=forward_molecule, 
        reverse_molecule=reverse_molecule,
        user_parameters=config.get("geom_opt", {}), 
        num_tasks=config.get("submit", {"ntasks":2})["ntasks"]
        )
    
    # save the 3 molecules (forward, backward, and TS) in a dedicated folder
    os.mkdir("./final_structures")
    os.chdir("./final_structures")
    forward_optimized.to("FORWARD.xyz")
    reverse_optimized.to("REVERSE.xyz")
    transition_state.to("TRANSITION_STATE.xyz")

    print("Program Finished Gracefully.\n Have a Nice Day :)")



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
    ts2rxn(config)