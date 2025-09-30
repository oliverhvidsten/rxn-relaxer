from rxnrlx.jaguar.read_files import get_mols_from_irc, get_energy_from_file
import os

DIR_PATH = os.path.dirname(__file__)


def test_get_mols__simple():
    """
    Test to see if the 
    """
    try:
        forward, backward = get_mols_from_irc(f"{DIR_PATH}/inputs/irc.out", 7)
    except:
        assert False
    
    assert (forward is not None) and (backward is not None)


def test_get_energy():
    """
    Ensure that correct energy is found in the file
    """

    energy = get_energy_from_file(f"{DIR_PATH}/inputs/energy_rev.out")

    assert energy == -799.720018


def test_verify_success():
    assert True