import numpy as np
from optparse import OptionParser
import math, random
#from mpl_toolkits.mplot3d import Axes3D
#import matplotlib.pyplot as plt

lattice_constants = {'Pt' : 0.39242}
lattice_constants['Au'] = 0.40782
lattice_constants['Fe'] = 0.28665
lattice_constants['Al'] = 0.40495
lattice_constants['Si'] = 0.54309
lattice_constants['Ti'] = 0.29508
lattice_constants['Pd'] = 0.38907
lattice_constants['Ag'] = 0.40853
lattice_constants['Cu'] = 0.36149

metal_radius = {'Pt' : 0.1385}
metal_radius['Au'] = 0.144
metal_radius['Fe'] = 0.126
metal_radius['Al'] = 0.143
metal_radius['Ti'] = 0.147
metal_radius['Pd'] = 0.137
metal_radius['Ag'] = 0.144
metal_radius['Cu'] = 0.128

parser=OptionParser()
parser.add_option("-o", "--output", action="store", type='string', dest="OutputFile", default='NP1.pdb', help="Name of the output file")
parser.add_option("-m", "--metal", action="store", type='string', dest="Metal", default='Au', help="Chemical symbol of the core")
parser.add_option("-r", "--radius", action="store", type='float', dest="Radius", default='3.0', help="Radius of the spherical core (nm)")
(options, args) = parser.parse_args()
outname_opt = options.OutputFile
metal_opt = options.Metal
radius_opt = options.Radius
fcc = False
hollow = True
solid = False

def center(objeto):
    COM = np.average(objeto, axis=0)
    for i in range(len(objeto)):
        objeto[i,:] = objeto[i,:] - COM
    return objeto

def fcc_solid_sphere():
    const = lattice_constants[metal_opt]
    cells_per_side = int((((2*radius_opt)//const)+1)//2*2+1)
    N_unit_cells = cells_per_side**3
    N_beads = N_unit_cells * 14
    fcc_block = np.array([])

    for i in range(cells_per_side):
        for j in range(cells_per_side):
            for k in range(cells_per_side):

                fcc_block = np.append(fcc_block, [i,j,k,i+1,j,k,i,j+1,k,i,j,k+1,i+1,j+1,k,i+1,j,k+1,i,j+1,k+1,i+1,j+1,k+1])
                fcc_block = np.append(fcc_block, [i,j+0.5,k+0.5,i+0.5,j,k+0.5,i+0.5,j+0.5,k,i+1,j+0.5,k+0.5,i+0.5,j+1,k+0.5,i+0.5,j+0.5,k+1])

    fcc_block = fcc_block * const
    fcc_block = fcc_block.reshape((N_beads,3))
    fcc_block = np.unique(fcc_block, axis=0)
    fcc_block = center(fcc_block)

    fcc_sphere=fcc_block[np.linalg.norm(fcc_block, axis=1)<= radius_opt]

    return fcc_sphere

def hollow_sphere():
    A = 4 * radius_opt**2
    a = metal_radius[metal_opt]**2
    N = int(A//a)
    return fibonacci_sphere(N)

def fibonacci_sphere(samples):
    rnd = 1.

    points = []
    offset = 2./samples
    increment = math.pi * (3. - math.sqrt(5.));

    for i in range(samples):
        y = ((i * offset) - 1) + (offset / 2);
        r = math.sqrt(1 - pow(y,2))

        phi = ((i + rnd) % samples) * increment

        x = math.cos(phi) * r
        z = math.sin(phi) * r

        points.append([x,y,z])

    return radius_opt*np.array(points)

def solid_sphere():

def print_xyz(coords):
    coords = coords * 10
    print(str(len(coords)))
    print(" ")
    for i in range(len(coords)):
        print(metal_opt + ' {:.3f}  {:.3f}  {:.3f}'. format(coords[i,0], coords[i,1], coords[i,2]))

if fcc == True:
    print_xyz(fcc_solid_sphere())
elif hollow == True:
    print_xyz(hollow_sphere())
elif solid == True:
    print_xyz(solid_sphere)
