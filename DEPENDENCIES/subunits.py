import numpy as np

class Block:
    def __init__(self, ndx_S=0, ndx_Au=0, ndx_C=0, ndx_H=0):
        self.S = ndx_S
        self.Au = ndx_Au
        self.C = ndx_C
        self.H = ndx_H

class Staple:
    def __init__(self, tipo='UNK', ndx_S=0, ndx_Au=0, ndx_Au_l=0, ndx_C=0, ndx_H=0):
        self.tipo = tipo
        self.S = np.array(ndx_S, dtype='int')
        self.Au_l = np.array(ndx_Au_l, dtype='int')
        self.Au = np.array(ndx_Au, dtype='int')
        self.Au_s = np.array(list(set(list(self.Au))-set(list(self.Au_l))), dtype='int')
        self.C = np.array(ndx_C, dtype='int')
        self.H = np.array(ndx_H, dtype='int')

    def change_tipo(self, new_tipo):
        self.tipo = new_tipo

class Residue:
    def __init__(self, restype='UNK', ndx=[]):
        self.restype = restype
        self.ndx = ndx
