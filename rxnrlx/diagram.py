"""
TODO: 
This file should take in molecules and energies in a form 
similar to what refine.py outputs and output a reaction diagram.

It should be possible to chain reactions together as well
"""
import os, sys, yaml
from pymatgen.core.structure import Molecule
from typing import Union

import matplotlib.pyplot as plt
from energydiagram import ED


def create_diagram(config:dict):
    """
    --- Example Config File ---
    info:
        subfolders: [subfolder1, subfolder2]
        order:                          # this specifies which way each reaction runs 
            subfolder1: reverse
            subfolder2: forward
        order: exergonic          # a different way to specify which wa
        omit_structures:
            subfolder1: [FORWARD.xyz]
            subfolder2: []
        use_refined_structures: True    # default is True (if available)
        energy_unit: eV           # accepted options: [eV, kcal]
    """

    """
    Take different relaxation jobs and put them together in one big pathway
    The exact way they they are spliced together is to be speicified in the config file
    
    At the end, write a reaction diagram PNG and save the structures with their new names and ordering
    """

    # Get list of dictionaries from specified information
    structure_list = prepare_path(config.get("info", {}))


    # Create a new directory to save this information in
    os.makedirs(f"./full_path")
    os.chdir(f"./full_path")

    # Create energy plot and save it to a .png
    draw_diagram(structure_list)

    # Save all of the molecules with their new names
    os.makedirs(f"./structures")
    os.chdir(f"./structures")
    
    for mol_dict in structure_list:
        mol_dict["molecule"].to(f"{mol_dict['name']}.xyz")


def draw_diagram(structure_list):

    diagram = ED()
    #diagram.dimension = 100
    plt.rcParams["font.family"] = "serif"
    plt.rcParams["font.family"] = ["Georgia"]

    for i, mol_dict in enumerate(structure_list):
        diagram.add_level(
            energy=mol_dict["energy"],
            bottom_text=mol_dict["name"]
        )
        if i != 0:
            diagram.add_link(i-1, i)
    
    diagram.plot(ylabel="Energy [$eV$]")

    diagram.dimension = 0.5

    diagram.plot(ylabel="Energy [$eV$]")
    #diagram.fig.set_figheight(10)
    #plt.show()
    plt.savefig(f"./reaction_diagram.png", dpi=300, bbox_inches='tight')



def prepare_path(info:dict) -> list[dict]:

    # initialize full path list to be added to
    full_path = list()

    ts_counter = 0
    stable_counter = 0

    # open the subfolders in order
    for folder in info.get("subfolder", {}):

        # get the path to the structures of interest
        if os.path.exists(f"./{folder}/refine_structures/final_structures") and info.get("use_refined_structures", True):
            structure_dir = f"./{folder}/refine_structures/final_structures"
            backup_structure_dir = f"./{folder}/final_structures"
        else:
            structure_dir = f"./{folder}/final_structures"
            backup_structure_dir = None

        # get energy dict
        with open(f"./{folder}/refined_structures/energy.yaml", "r") as f:
            energy_dict = yaml.safe_load(f)

        # change energy values to requested 
        energy_dict = convert_energy_units(energy_dict, info.get("energy_unit", "eV"))



        structure_order = get_order(
            order_parameter=info.get("order"), 
            omit_parameter=info.get("omit_structures"), 
            energy_dict=energy_dict,
            folder=folder
        )

        for filename in structure_order:
            mol_dict = prepare_structure_dict(
                structure_dir=structure_dir, 
                backup_structure_dir=backup_structure_dir,
                energy_dict=energy_dict,
                filename=filename
            )
            mol_dict["name"], ts_counter, stable_counter = get_name(mol_dict, ts_counter, stable_counter)
            
            full_path.append(mol_dict)
            
    return full_path


def convert_energy_units(energy_dict, units):
    """
    Convert energy values from Hartrees to a more human interpretable unit 
    """

    if units == "eV": 
        conversion = 27.2114
    elif units == "kcal":
        conversion = 627.5095
    else:
        raise Exception("Unrecognized Energy Unit: Please choose between: 'eV' and 'kcal'")
    
    for mol in ["forward", "reverse", "transition_state"]:
        energy_dict[mol] = energy_dict[mol] * conversion

    return energy_dict


def get_name(mol_dict:dict, ts_counter:int, stable_counter:int):
    """
    Creates the name of the structure
    """
    if mol_dict["type"] == "ts":
        ts_counter += 1
        name = f"TS{ts_counter}"
    elif mol_dict["type"] == "stable":
        name = f"M{stable_counter}"
        stable_counter += 1

    return name, ts_counter, stable_counter

def get_order(order_parameter, omit_parameter, energy_dict, folder) -> list:
    """
    Get the order in which the molecules will be added to the reaction flow diagram
    Remove any molecules that have been specified to be omitted
    """

    if isinstance(order_parameter, dict):
        # if it is a dictionary, the user has specified the order of each reaction step
        direction = order_parameter.get(folder)
        
    elif isinstance(order_parameter, str):
        # if it is a string, the user has specific if they want all of the reactions to flow exergonically or endergonically
        if order_parameter == "exergonic":
            if energy_dict["reverse"] > energy_dict["forward"]:
                 direction = "forward"
            else:
                 direction = "reverse"

        elif order_parameter == "endergonic":
            if energy_dict["forward"] > energy_dict["reverse"]:
                 direction = "forward"
            else:
                 direction = "reverse"
        else: 
            raise Exception(f"Unrecognized Request: '{order_parameter}' is not a valid option, please select 'exergonic' or 'endergonic'.")

    
    if direction == "forward":
        structure_order = ["REVERSE.xyz", "TRANSITION_STATE.xyz", "FORWARD.xyz"]
    elif direction == "reverse":
        structure_order = ["FORWARD.xyz", "TRANSITION_STATE.xyz", "REVERSE.xyz"]
    else:
        raise Exception(f"Unrecognized Direction: '{direction}' is not a valid option, please select 'forward' or 'reverse'.")
    

    # check if any of the files should be omitted
    omit_files = omit_parameter.get(folder, [])
    for file in omit_files:
        if file in structure_order:
            structure_order.remove(file)
        else:
            raise Exception(f"Unrecognized File: '{file}' is not contained with in the valid structures {structure_order}")

    # Return final structure order       
    return structure_order


def prepare_structure_dict(structure_dir:str, backup_structure_dir:Union[str, None], energy_dict:dict, filename:str) -> dict:
    """
    Create the dictionary that will be passed into the structures list
    """
    # Check if the file exists with the main path
    if os.path.exists(f"{structure_dir}/{filename}"):
        mol = Molecule.from_file(f"{structure_dir}/{filename}")
    elif (backup_structure_dir is not None) and os.path.exists(f"{backup_structure_dir}/{filename}"):
        print(f"\nWARNING: Defaulting to pre-refined structure because file: '{structure_dir}/{filename}' does not exist.\n")
        mol = Molecule.from_file(f"{backup_structure_dir}/{filename}")
    else:
        raise Exception(f"Have all structures been optimized? '{structure_dir}/{filename}' does not exist.")

    molecule_name, _ = filename.split(".")

    # get type of structure (stable geometry or transition state)
    if molecule_name == "forward" or molecule_name == "reverse":
        molecule_type = "stable"
    elif molecule_name == "transition_state":
        molecule_type = "ts"
    else:
        raise Exception("Unrecongized Molecule Name: This is an unexpected error, please ensure you have not changed any filenames.")

    energy = energy_dict[molecule_name]



    return {
        "molecule": mol,
        "energy": energy, # convert to 
        "type": molecule_type
    }
    




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
    create_diagram(config)