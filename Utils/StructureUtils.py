from Bio.PDB import NeighborSearch

def get_neighbors(structure, radius=20):
    all_neighbors = _get_all_neighbors(structure, radius)
    return all_neighbors
    #neighbors_from_different_chains = _get_neighbors_from_different_chains(all_neighbors)
    #return _get_neighbors_by_seq_col(neighbors_from_different_chains)

def _get_all_neighbors(structure, radius):
    atoms = [atom for atom in structure.get_atoms()]
    neighbors_search = NeighborSearch(atoms)
    return neighbors_search.search_all(radius, level='R')

def _get_neighbors_from_different_chains(neighbors):
    None
