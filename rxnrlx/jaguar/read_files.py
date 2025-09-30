from pymatgen.core.structure import Molecule

import re

def get_energy_from_file(outfile:str) -> float:
    """ Get value of gibbs energy from energy output file """

    with open(outfile, "r") as f:
        lines = f.readlines()

    for line in lines:
        if "Total Gibbs free energy" in line:
            _, energy = line.split(":")
            energy, _ = energy.split()
    
    return float(energy)




def get_mols_from_irc(outfile:str, num_atoms:int) -> tuple[Molecule, Molecule]:
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


def get_mol_from_opt(outfile:str, num_atoms:int) -> Molecule:
    """ Get Molecule out of a optimizaiton job (TS or Stable Geometry)"""
    
    with open(outfile, "r") as f:
        lines = f.readlines()

    return find_molecule_in_section(lines, len(lines)-1, num_atoms)



def find_molecule_in_section(lines, starting_place, num_atoms) -> Molecule:
    """ Find the first relaxed molecule definition to appear before the given line index """

    # Find the geometry header line
    found = False
    i = starting_place
    while not found:
        pattern = re.compile(r"atom\s+x\s+y\s+z")
        if re.search(pattern, lines[i]):
            found = True
        i -= 1 # iterating up the file now
    
    i+=2 # i should now hold the first line of the atoms
    
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


def verify_success(outfile, name):
    """
    Check the outfile for language that verifies that the job was completed successfully
    """

    with open(outfile, "r") as f:
        lines = f.readlines()

    success_pattern = re.compile(rf'Job {name} completed on')

    return re.match(success_pattern, lines[-1])