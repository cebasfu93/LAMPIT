import numpy as np
from optparse import OptionParser
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt

lattice_constants = {'Pt' : 0.39242}
lattice_constants['Au'] = 0.40782
lattice_constants['Fe'] = 0.28665
lattice_constants['Al'] = 0.40495
lattice_constants['Si'] = 0.54309
lattice_constants['Ti'] = 0.29508
lattice_constants['Pd'] = 0.38907
lattice_constants['Ag'] = 0.40853
lattice_constants['Cu'] = 0.36149

parser=OptionParser()
parser.add_option("-o", "--output", action="store", type='string', dest="OutputFile", default='NP1.pdb', help="Name of the output file")
parser.add_option("-m", "--metal", action="store", type='string', dest="Metal", default='Au', help="Chemical symbol of the core")
parser.add_option("-r", "--radius", action="store", type='float', dest="Radius", default='3.0', help="Radius of the spherical core (nm)")
(options, args) = parser.parse_args()
outname_opt = options.OutputFile
metal_opt = options.Metal
radius_opt = options.Radius
const = lattice_constants[metal_opt]

def center(objeto):
    COM = np.average(objeto, axis=0)
    for i in range(len(objeto)):
        objeto[i,:] = objeto[i,:] - COM
    return objeto

cells_per_side = int((((2*radius_opt)//const)+1)//2*2+1)
N_unit_cells = cells_per_side**3
N_beads = N_unit_cells * 14
metal_block = np.array([])

for i in range(cells_per_side):
    for j in range(cells_per_side):
        for k in range(cells_per_side):

            metal_block = np.append(metal_block, [i,j,k,i+1,j,k,i,j+1,k,i,j,k+1,i+1,j+1,k,i+1,j,k+1,i,j+1,k+1,i+1,j+1,k+1])
            metal_block = np.append(metal_block, [i,j+0.5,k+0.5,i+0.5,j,k+0.5,i+0.5,j+0.5,k,i+1,j+0.5,k+0.5,i+0.5,j+1,k+0.5,i+0.5,j+0.5,k+1])

metal_block = metal_block * const
metal_block = metal_block.reshape((N_beads,3))
metal_block = np.unique(metal_block, axis=0)
metal_block = center(metal_block)

metal_sphere=metal_block[np.linalg.norm(metal_block, axis=1)<= radius_opt]

metal_sphere = metal_sphere * 10
print(str(len(metal_sphere)))
print(" ")
for i in range(len(metal_sphere)):
    print('Pt  {:.3f}  {:.3f}  {:.3f}'. format(metal_sphere[i,0], metal_sphere[i,1], metal_sphere[i,2]))
