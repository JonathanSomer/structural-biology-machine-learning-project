from Bio.PDB import PDBParser, PPBuilder
from Bio.PDB.Residue import Residue
from Bio.PDB.Structure import Structure
from Bio.PDB.Superimposer import Superimposer
from Bio.PDB.PDBExceptions import PDBConstructionException
from Bio import pairwise2

pdb_parser = PDBParser()


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


def get_rmsd_of_super_imposed_complex(complex1, complex2):
    super_imposer = Superimposer()
    complex1_atom_list, complex2_atom_list = align_complexes_atoms_lsts(complex1, complex2)
    super_imposer.set_atoms(complex1_atom_list, complex2_atom_list)
    return super_imposer.rms


def align_complexes_atoms_lsts(complex1, complex2):
    """
    Gets two complexes complex1, complex2
    Return two lists containing only the pairs of atoms that appear in complex1 and complex2
    """
    complex1_atoms_aligned = []
    complex2_atoms_aligned = []
    for residue1, residue2 in zip(complex1.get_residues(),
                                  complex2.get_residues()):
        residue1_atoms_aligned, residue2_atoms_aligned = _align_residues_atom_lsts(residue1, residue2)
        complex1_atoms_aligned += residue1_atoms_aligned
        complex2_atoms_aligned += residue2_atoms_aligned
    return complex1_atoms_aligned, complex2_atoms_aligned


def _align_residues_atom_lsts(residue1, residue2):
    """
    Gets two residues residue1, residue2
    Return two lists containing only the pairs of atoms that appear in residue1 and residue2
    """
    residue1_atoms_aligned = []
    residue2_atoms_aligned = []
    if residue1.resname != residue2.resname:
        raise PDBConstructionException("""expected residues to be of same AA. 
                                       Got {} on {} and {} on {}""".format(residue1.resname, residue1.id[1],
                                                                           residue2.resname, residue2.id[1]))
    for a1 in residue1:
        for a2 in residue2:
            if a1.fullname == a2.fullname:
                residue1_atoms_aligned.append(a1)
                residue2_atoms_aligned.append(a2)
    return residue1_atoms_aligned, residue2_atoms_aligned

def get_joint_positions_from_pdbs(path_to_pdb1, path_to_pdb2):
    """get paths to two pdb files
    return a list oa list of tuples of the identical residues.
    """
    alignments = pairwise2.align.globalxx(
        get_structure_sequence(pdb_parser.get_structure("1", path_to_pdb1)),
        get_structure_sequence(pdb_parser.get_structure("2", path_to_pdb2)))
    return get_joint_positions_from_alignment(*alignments[0])


def get_joint_positions_from_alignment(align1, align2, score, begin, end):
    """format_alignment(align1, align2, score, begin, end) -> string
    Format the alignment prettily into a list of tuples of the identical residues.
    """
    lst = []
    i, align1_index, align2_index = 0, 0, 0
    while i < len(align1):
        if align1[i] == align2[i]:
            lst.append((align1_index, align2_index))
        if align1[i] != "-":
            align1_index += 1
        if align2[i] != "-":
            align2_index += 1
        i += 1
    return lst