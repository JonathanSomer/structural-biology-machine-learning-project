from Bio.PDB import PPBuilder

def get_sequence(struct):
    ppb = PPBuilder()
    return ppb.build_peptides(struct)[0].get_sequence()