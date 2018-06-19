from Bio.PDB import PPBuilder
from Bio.PDB.Residue import Residue
from Bio.PDB.Structure import Structure


def is_protein_residue(residue):
    # type: (Residue) -> bool
    """
    Checks if the residue comes from a protein
    :param residue: the given residue
    :return: is from protein
    """
    return residue.id[0] == ' '


def get_residue_global_id(residue):
    # type: (Residue) -> int
    """
    Gets the global index in the structure
    returns -1 if residue is not s protein residue
    :param residue: Bio.PDB.Residue.Residue object
    :return: global id (integer)
    """
    if not is_protein_residue(residue):
        return -1
    struct = residue.get_parent().get_parent().get_parent()
    global_id = 0
    for r in struct.get_residues():
        if not is_protein_residue(r):
            continue
        # stop when we get to the residue in the structure
        if r.full_id == residue.full_id:
            break
        global_id += 1
    return global_id


def get_structure_sequence(struct):
    # type: (Structure) -> str
    """
    Gets the structure sequence using PPBuilder
    :param struct: Structure object
    :return: struct sequence
    """
    ppb = PPBuilder()
    return ''.join([str(pp.get_sequence()) for pp in ppb.build_peptides(struct)])
