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

def jaguar_input(file_name:str, structure:Molecule, parameters:dict={}):
    """
    Create an input file for a simple jaguar job for one structure
    """
    # Create gen section
    gen_section = create_gen_section(parameters)

    # Create zmat section
    zmat_section = structure.get_zmatrix()

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

    

    """
    gen_section = [
        "&gen",
        f"igeopt = 2",
        f"inhess = {nondefault_parameters.get('inhess', 4)}",
        f"epsout = 18.5",
        f"isolv = 7",
        f"maxitg = 300",
        f"isymm = 0",
        f"basis = DEF2-SVPD",
        f"maxit = 300",
        f"ip472 = 2",
        f"dftname = wB97X-D",
        f"nogas = 2",
        f"ip175 = 2",
        f"iacc = 2",
        f"&",
    ]
    &gen
    """

