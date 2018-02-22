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
parser.add_option("-m", "--metal", action="store", type='string', dest="Metal", default='Pt', help="Chemical symbol of the core")
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

cells_per_side=int((((2*radius_opt)//const)+1)//2*2+1)
metal_block = np.array([0.0, 0.0, 0.0])
for i in range(cells_per_side):
    for j in range(cells_per_side):
        for k in range(cells_per_side):
            xyz1 = np.array([i*const, j*const, k*const])
            xyz2 = np.array([(i+1)*const, j*const, k*const])
            xyz3 = np.array([i*const, (j+1)*const, k*const])
            xyz4 = np.array([i*const, j*const, (k+1)*const])
            xyz5 = np.array([(i+1)*const, (j+1)*const, k*const])
            xyz6 = np.array([(i+1)*const, j*const, (k+1)*const])
            xyz7 = np.array([i*const, (j+1)*const, (k+1)*const])
            xyz8 = np.array([(i+1)*const, (j+1)*const, (k+1)*const])

            xyz9 = np.array([i*const, (j+0.5)*const, (k+0.5)*const])
            xyz10 = np.array([(i+0.5)*const, j*const, (k+0.5)*const])
            xyz11 = np.array([(i+0.5)*const, (j+0.5)*const, k*const])
            xyz12 = np.array([(i+1)*const, (j+0.5)*const, (k+0.5)*const])
            xyz13 = np.array([(i+0.5)*const, (j+1)*const, (k+0.5)*const])
            xyz14 = np.array([(i+0.5)*const, (j+0.5)*const, (k+1)*const])

            metal_block = np.vstack((metal_block, xyz1))
            metal_block = np.vstack((metal_block, xyz2))
            metal_block = np.vstack((metal_block, xyz3))
            metal_block = np.vstack((metal_block, xyz4))
            metal_block = np.vstack((metal_block, xyz5))
            metal_block = np.vstack((metal_block, xyz6))
            metal_block = np.vstack((metal_block, xyz7))
            metal_block = np.vstack((metal_block, xyz8))
            metal_block = np.vstack((metal_block, xyz9))
            metal_block = np.vstack((metal_block, xyz10))
            metal_block = np.vstack((metal_block, xyz11))
            metal_block = np.vstack((metal_block, xyz12))
            metal_block = np.vstack((metal_block, xyz13))
            metal_block = np.vstack((metal_block, xyz14))
metal_block = np.unique(metal_block, axis=0)
metal_block = center(metal_block)

"""
metal_block = metal_block * 10
print(str(len(metal_block)))
print(" ")
for i in range(len(metal_block)):
    print("Pt  {:.3f}  {:.3f}  {:.3f}". format(metal_block[i,0], metal_block[i,1], metal_block[i,2]))
"""
metal_sphere=np.array([0.0, 0.0, 0.0])
for i in range(len(metal_block[:,0])):
    if np.linalg.norm(metal_block[i,:]) <= radius_opt:
        metal_sphere = np.vstack((metal_sphere, metal_block[i,:]))
metal_sphere = np.delete(metal_sphere, 0, 0)

metal_sphere = metal_sphere * 10
print(str(len(metal_sphere)))
print(" ")
for i in range(len(metal_sphere)):
    print('Pt  {:.3f}  {:.3f}  {:.3f}'. format(metal_sphere[i,0], metal_sphere[i,1], metal_sphere[i,2]))
