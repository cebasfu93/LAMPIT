class Residue:
    def __init__(self, number, type, gold_atoms, sulphur_atoms):
        self.restype=type
        self.resnum=number
        self.au_atoms=gold_atoms
        self.s_atoms=sulphur_atoms
        self.N_atoms=len(self.s_atoms)+len(self.au_atoms)

    def change_type(self, new_type):
        self.restype=new_type
