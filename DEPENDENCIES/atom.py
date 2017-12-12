class Residue:
    def __init__(self, number, type, gold_atoms, sulphur_atoms, gold_n_bonds):
        self.restype=type
        self.resnum=number
        self.au_atoms=gold_atoms
        self.s_atoms=sulphur_atoms
        self.N_atoms=len(self.s_atoms)+len(self.au_atoms)
        self.au_n_bonds=gold_n_bonds

    def change_type(self, new_type):
        self.restype=new_type
