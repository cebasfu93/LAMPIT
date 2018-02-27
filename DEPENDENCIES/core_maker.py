import numpy as np
from optparse import OptionParser
import math, random
from  scipy.spatial.distance import cdist
#from mpl_toolkits.mplot3d import Axes3D
#import matplotlib.pyplot as plt

lattice_constants = {'Pt' : 0.39242} #in nm
lattice_constants['Au'] = 0.40782
lattice_constants['Fe'] = 0.28665
lattice_constants['Al'] = 0.40495
lattice_constants['Si'] = 0.54309
lattice_constants['Ti'] = 0.29508
lattice_constants['Pd'] = 0.38907
lattice_constants['Ag'] = 0.40853
lattice_constants['Cu'] = 0.36149

metal_radius = {'Pt' : 0.1385} #in nm
metal_radius['Au'] = 0.144
metal_radius['Fe'] = 0.126
metal_radius['Al'] = 0.143
metal_radius['Ti'] = 0.147
metal_radius['Pd'] = 0.137
metal_radius['Ag'] = 0.144
metal_radius['Cu'] = 0.128

metal_radius['S'] = 0.102

parser=OptionParser()
parser.add_option("-o", "--output", action="store", type='string', dest="OutputFile", default='NP1.pdb', help="Name of the output file")
parser.add_option("-m", "--metal", action="store", type='string', dest="Metal", default='Au', help="Chemical symbol of the core")
parser.add_option("-n", "--ligands", action="store", type='int', dest="NumberLigands", default='60', help="Number of ligands to place around the sphere")
parser.add_option("-r", "--radius", action="store", type='float', dest="Radius", default='3.0', help="Radius of the spherical core (nm)")
parser.add_option("-t", "--type", action="store", type='string', dest="SphereType", default='hollow', help="Methods used to build the sphere. Valid options: fcc, hollow, solid, semisolid. When solid is selected, the output radius is not the same as the requested")
parser.add_option("-d", "--thickness", action="store", type='float', dest="Thickness", default='1.4', help="Thickness of the shell (nm) when using option semisolid")
(options, args) = parser.parse_args()
outname_opt = options.OutputFile
metal_opt = options.Metal
lignum_opt = options.NumberLigands
radius_opt = options.Radius
sphere_opt = options.SphereType
thickness_opt = options.Thickness

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

    fcc_sphere = np.vstack((fcc_sphere, put_staples(fcc_sphere, radius_opt)))

    return fcc_sphere

def hollow_sphere(radius, samples):
    return (radius-metal_radius[metal_opt])*fibonacci_sphere(samples)

def get_N(radius):
    A = 4 * radius**2
    a = metal_radius[metal_opt]**2
    N = int(A//a)
    return N

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

    return np.array(points)

def solid_sphere():
    D = 4 / (lattice_constants[metal_opt]**3) #atoms/nm3

    N_shells = int(round(radius_opt / (2*metal_radius[metal_opt])))+1

    real_radius = (2*N_shells-1)*metal_radius[metal_opt]
    V = 4.0 * math.pi * real_radius**3 /3.0 #nm3

    N_atoms = D * V
    atoms_shells = []
    for i in range(N_shells):
        atoms_shells.append(get_N(2*i*metal_radius[metal_opt]))

    atoms_shells[0] = 1
    delta = N_atoms - np.sum(np.array(atoms_shells))
    points = np.array([0.0, 0.0, 0.0])

    for i in range(N_shells):
        atoms_shells[i] += int(3*delta*(i**2)/N_shells**3)
        shell = hollow_sphere(2*i*metal_radius[metal_opt], atoms_shells[i])
        if i == (N_shells-1):
            S_atoms = put_staples(shell, real_radius-metal_radius[metal_opt])
        points = np.vstack((points, shell))
    points = np.vstack((points, S_atoms))
    points = np.delete(points, 0, 0)

    print("The output radius is {:.3f} nm".format(real_radius))
    print("The theoretical density is {:.3f} atoms/nm^3".format(D))
    print("The achieved density is {:.3f} atoms/nm^3".format(len(points)/V))

    return points

def semi_solid_sphere(thick):
    if thick >= radius_opt:
        print("Thickness bigger than NP radius")

    D = 4 / (lattice_constants[metal_opt]**3) #atoms/nm3

    N_shells = int(round(thick / (2*metal_radius[metal_opt])))+1

    real_shell = (2*N_shells)*metal_radius[metal_opt]
    V = 4.0 * math.pi / 3.0 *(radius_opt**3 - real_shell**3)

    N_atoms = D * V
    atoms_shells = []
    for i in range(N_shells):
        atoms_shells.append(get_N(radius_opt - real_shell + (2*i+1)*metal_radius[metal_opt]))

    delta = N_atoms - np.sum(np.array(atoms_shells))
    points = np.array([0.0, 0.0, 0.0])

    for i in range(N_shells):
        atoms_shells[i] += int(3*delta*((i+0.5)**2)/(N_shells**3))
        shell = hollow_sphere(radius_opt - real_shell + (2*i+1)*metal_radius[metal_opt], atoms_shells[i])
        if i == (N_shells-1):
            S_atoms = put_staples(shell, radius_opt-metal_radius[metal_opt])
        points = np.vstack((points, shell))
    points = np.vstack((points, S_atoms))
    points = np.delete(points, 0, 0)

    print("The theoretical density is {:.3f} atoms/nm^3".format(D))
    print("The achieved density is {:.3f} atoms/nm^3".format(len(points)/V))
    return points

def put_staples(shell, radius):
    S_atoms = hollow_sphere(radius, lignum_opt)
    distances = cdist(S_atoms, shell)
    mins = np.argmin(distances, axis=1)
    for i in range(len(S_atoms)):
        norma=np.linalg.norm(shell[mins[i]])
        scaling = (norma + (metal_radius[metal_opt]+metal_radius['S'])/2)/norma
        S_atoms[i,:] = scaling*shell[mins[i]]
    return S_atoms

def print_xyz(coords):
    coords = coords * 10
    output = open(outname_opt, "a")
    output.write(str(len(coords)) + "\n\n")
    N_metal = len(coords)-lignum_opt
    for i in range(N_metal):
        output.write(metal_opt + '{:.3f}'.format(coords[i,0]).rjust(10) + "{:.3f}".format(coords[i,1]).rjust(10) + "{:.3f}".format(coords[i,2]).rjust(10) + "\n")
    for i in range(N_metal, len(coords)):
        output.write('S' + '{:.3f}'.format(coords[i,0]).rjust(10) + "{:.3f}".format(coords[i,1]).rjust(10) + "{:.3f}".format(coords[i,2]).rjust(10) + "\n")
    output.close()

if sphere_opt == "fcc":
    print_xyz(fcc_solid_sphere())
elif sphere_opt == "hollow":
    points = hollow_sphere(radius_opt, get_N(radius_opt))
    NP = np.vstack((points, put_staples(points, radius_opt-metal_radius[metal_opt])))
    print_xyz(NP)
elif sphere_opt == "solid":
    print_xyz(solid_sphere())
elif sphere_opt == "semisolid":
    print_xyz(semi_solid_sphere(thickness_opt))
