from pymatgen.core.structure import Molecule


def create_gen_section(parameters:dict={}) -> str:
    """
    Use the user-specified paramters to create jaguar's gen section
    """
    # Put specifications in 
    gen_section = ["&gen"]
    for key,value in parameters.items():
        gen_section.append(f"{key} = {value}")
    gen_section.append("&")

    return "\n".join(gen_section)


def create_zmat_section(structure:Molecule):
    """
    Use the specified structure to create the zmat section of the input file
    """
    zmat_section = ["&zmat"]
    for i, site in enumerate(structure.sites):
        x,y,z = site.coords
        symb = f"{site.species_string}{i}"
        zmat_section.append(f"{symb:<5s} {x: .9f} {y: .9f} {z: .9f}")
    zmat_section.append("&")

    return "\n".join(zmat_section)

def jaguar_input(file_name:str, structure:Molecule, parameters:dict={}):
    """
    Create an input file for a simple jaguar job for one structure
    """
    # Create gen section
    gen_section = create_gen_section(parameters)

    # Create zmat section
    zmat_section = create_zmat_section(structure)

    with open(file_name, "w") as f:
        f.write(f"{gen_section}\n{zmat_section}")
    

def multi_species_jaguar_input(structure:list[Molecule], parameters:dict={}):
    """
    Some jobs like QST transition state optimizations might need multiple zmat sections
    Currently, nothing in this workflow needs this so it remains unimplemnted for now.
    """
    # Create gen seciton
    gen_section = create_gen_section(parameters)
    # TODO: create multiple zmat sections (check implementation requirements)
    raise NotImplementedError

