import numpy as np
import matplotlib.pyplot as plt
from optparse import OptionParser
from sklearn.preprocessing import normalize

parser=OptionParser()
parser.add_option("-o", "--output", action="store", type='string', dest="OutputFile", default='Au144L1-60.xyz', help="Name of the output file")
parser.add_option("-c", "--core", action="store", type='string', dest="CoreFile", default='Au144_MD6341surface1.xyz', help="Name of the xyz core file (includes S)")
parser.add_option("-l", "--ligand", action="store", type='string', dest="LigandFile", default='ligand.mol2', help="Name of the mol2 uncapped ligand file")
parser.add_option("-r", "--resname", action="store", type='string', dest="ResidueName", default='UNK', help="3-letter name of the coating ligand")
(options, args)= parser.parse_args()
outname=options.OutputFile
corename=options.CoreFile
ligname=options.LigandFile
res_name=options.ResidueName

def init_lig_mol2(ligand_fname):
    lig_file=np.genfromtxt(ligand_fname, delimiter='\n', dtype='str')
    for i in range(len(lig_file)):
        if "@<TRIPOS>ATOM" in lig_file[i]:
            ini=i+1
        if "@<TRIPOS>BOND" in lig_file[i]:
            fin=i
            break
    lig_file=lig_file[ini:fin]
    n_at=len(lig_file)
    names=[]
    types=[]
    xyz_lig=np.zeros((n_at,3))
    for i in range(n_at):
        at_file=lig_file[i].split()
        names.append(at_file[1])
        types.append(at_file[5])
        xyz_lig[i,:]=at_file[2:5]
    return xyz_lig, names, types

def init_core_xyz(core_fname):
    core_file=np.genfromtxt(core_fname, skip_header=2, dtype='str')
    name_core=core_file[:,0]
    xyz_core=core_file[:,1:].astype('float')
    return xyz_core, name_core

core, names_core=init_core_xyz(corename)
lig, names_lig, types_lig=init_lig_mol2(ligname)

N_at_lig=len(lig[:,0])
N_S=len(names_core[names_core=='S'])
d_CS=1.83 #Angstroms (wikipedia)

def get_COM(xyz, objeto_COM):
    return np.average(xyz[names_core==objeto_COM,:], axis=0)
def get_new_C(C_of_mass, xyz_s, dist):
    COM_S=np.subtract(xyz_s, C_of_mass)
    mag_COM_S=np.linalg.norm(COM_S)
    const=((mag_COM_S+dist)/mag_COM_S)**0.5
    new_C=C_of_mass+const*COM_S
    return new_C
def get_rot_matrix(a, b):
    matrix=np.identity(3)
    v=np.cross(a,b)
    s=np.linalg.norm(v)
    c=np.dot(a,b)
    vx=np.array([[0, -v[2], v[1]],[v[2], 0, -v[0]],[-v[1], v[0], 0]])
    vx_2=np.dot(vx, vx)
    factor=(1-c)/s**2
    matrix=matrix+vx+vx_2*factor
    return normalize(matrix)
def place_ligands(C_of_mass, xyz_core, xyz_lig, terminal_ndx):
    N_at_lig_tot=N_at_lig*N_S
    xyz_S=xyz_core[np.where(names_core=='S')[0],:]
    new_ligs=np.zeros((N_at_lig_tot, 3))
    for i in range(N_S):
        old_C=lig[terminal_ndx,:]
        new_C=get_new_C(C_of_mass, xyz_S[i,:], d_CS)

        rot_matrix=get_rot_matrix(old_C, new_C)
        deltas=np.subtract(old_C,new_C)
        for j in range(N_at_lig):
            new_ligs[i*N_at_lig+j,:]=lig[j,:]+deltas
            new_ligs[i*N_at_lig+j,:]=np.dot(rot_matrix,lig[j,:])
    return new_ligs
def print_NP_pdb(name_lig, name_core, xyz_ligs, xyz_core, type_lig, filename):
    output=open(filename, "w")
    res=1
    at=0
    for i in range(len(xyz_core[:,0])):
        at+=1
        write_pdb_block(name_core[i], 'COR', xyz_core[i,:], res, at, filename)
    for i in range(N_S):
        res+=1
        for j in range(N_at_lig):
            at+=1
            write_pdb_block(name_lig[j], res_name, xyz_ligs[i*N_S+j,:], res, at, filename)
    output.close()
def write_pdb_block(atname, resname, xyz, resnum, atnum, file_name):
    xyz=np.round(xyz, decimals=4)
    coords=open(file_name, 'a')
    coords.write('ATOM'.ljust(6))
    coords.write(str(atnum).rjust(5))
    coords.write('  '+str(atname).ljust(3))
    coords.write(' '+str(resname).ljust(3))
    coords.write('  '+str(resnum).rjust(4))
    coords.write('   '+ str(xyz[0]).rjust(8))
    coords.write(str(xyz[1]).rjust(8))
    coords.write(str(xyz[2]).rjust(8)+"\n")
    coords.close()

COM=get_COM(core, 'Au')
ligands=place_ligands(COM, core, lig, 0)
print_NP_pdb(names_lig, names_core, ligands, core, types_lig, outname)
