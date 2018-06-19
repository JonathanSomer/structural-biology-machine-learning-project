from Bio.PDB.Residue import Residue


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
    :param residue: Bio.PDB.Residue.Residue object
    :return: global id (integer)
    """
    struct = residue.get_parent().get_parent().get_parent()
    global_id = residue.id[1]
    for chain in struct.get_chains():
        # break loop when we get to the residue's chain
        if chain.id == residue.full_id[2]:
            break

        global_id += len([r for r in chain.get_residues() if is_protein_residue(r)])
    return global_id
