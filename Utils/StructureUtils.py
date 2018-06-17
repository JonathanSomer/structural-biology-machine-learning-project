from Bio.PDB import NeighborSearch

def get_neighbors(structure, radius=20):
    atoms = [atom for atom in structure.get_atoms()]
    neighbors_search = NeighborSearch(atoms)
    all_neighbors = neighbors_search.search_all(radius, level='R')

