import numpy as np
import matplotlib.pyplot as plt
from optparse import OptionParser
import os
import sys
from  transformations import *

#Declaration of the flags to run the script
parser=OptionParser()
parser.add_option("-o", "--output", action="store", type='string', dest="OutputFile", default='NP1.xyz', help="Name of the output file")
parser.add_option("-c", "--core", action="store", type='string', dest="CoreFile", default='core1.pdb', help="Name of the pdb core file (includes S and anchor atom)")
parser.add_option("-l", "--ligand", action="store", type='string', dest="LigandFile", default='ligand.mol2', help="Name of the mol2 uncapped ligand file")
parser.add_option("-r", "--resname", action="store", type='string', dest="ResidueName", default='UNK', help="3-letter name to give to the coating ligand")
parser.add_option("-s", "--stones", action="append", dest="StonesIndexes", default=[], help="Indexes (starting from 0) of the atoms that will be aligned to the COM-S vector")
parser.add_option("-n", "--nucleus", action="store", type='string', dest="CoreAtomName", default='Au', help="Name of the atoms making up the core")
parser.add_option("-t", "--staple", action="store", type='string', dest="StapleAtomName", default='S', help="Name of the atoms making up the staples")
parser.add_option("-a", "--anchor", action="store", type='string', dest="AchorAtomName", default='C', help="Name of the anchor atom")
(options, args)= parser.parse_args()
outname_opt=options.OutputFile
corename_opt=options.CoreFile
ligname_opt=options.LigandFile
res_name_opt=options.ResidueName
stones_ndx_opt=np.array(list(map(int, options.StonesIndexes[0].split(','))))
core_at_name_opt=options.CoreAtomName
staple_at_name_opt=options.StapleAtomName
name_anchor_opt=options.AchorAtomName

def check_options(core_fname, ligand_fname, stones_fndx):
    #Checks that the required options exist and are consistent
    if not (os.path.isfile(core_fname)):
        sys.exit("Core file could not be found")
    if not (os.path.isfile(ligand_fname)):
        sys.exit("Ligand file could not be found")
    if not stones_fndx[-1]:
        sys.exit("No list of indexes was provided to do the roto-translation")

def init_lig_mol2(ligand_fname):
    #Imports ligand mol2 file. Returns xyz coordinates and names
    lig_file=np.genfromtxt(ligand_fname, delimiter='\n', dtype='str')
    for i in range(len(lig_file)):
        if "@<TRIPOS>ATOM" in lig_file[i]:
            ini=i+1
        if "@<TRIPOS>BOND" in lig_file[i]:
            fin=i
            break
    lig_file=lig_file[ini:fin]
    n_at=len(lig_file)
    names_lig_func=[]
    xyz_lig_func=np.zeros((n_at,3))
    for i in range(n_at):
        at_file=lig_file[i].split()
        names_lig_func.append(at_file[1])
        xyz_lig_func[i,:]=at_file[2:5]
    return np.array(xyz_lig_func), np.array(names_lig_func)

def init_core_pdb(core_fname):
    #Imports core pdb file. Returns xyz coordinates and names
    core_file=np.genfromtxt(core_fname, delimiter='\n', dtype=str)
    first=1
    N_core_func=0
    for i in range(len(core_file)):
        if "ATOM" in core_file[i] or "HETATM" in core_file[i]:
            N_core_func+=1
            if first==1:
                ini=i
                first=0
    core_file=core_file[ini:ini+N_core_func]
    names_core_func=[]
    xyz_core_func=np.zeros((N_core_func,3))
    for i in range(N_core_func):
        at_file=core_file[i].split()
        names_core_func.append(at_file[2])
        xyz_core_func[i,:]=at_file[5:8]
    xyz_core_func=np.array(xyz_core_func)
    names_core_func=np.array(names_core_func)

    COM=np.average(xyz_core_func[names_core_func=='Au',:], axis=0)

    for i in range(len(xyz_core_func[:,0])):
        xyz_core_func[i,:]=xyz_core_func[i,:]-COM
    return xyz_core_func, names_core_func

check_options(corename_opt, ligname_opt, stones_ndx_opt)
xyz_lig, names_lig=init_lig_mol2(ligname_opt)
xyz_core, names_core=init_core_pdb(corename_opt)

#Define number of atoms in a ligand, number of staples, and number of stones specified
N_at_lig=len(xyz_lig[:,0])
N_S=len(names_core[names_core==staple_at_name_opt])
N_stones=len(stones_ndx_opt)

def change_names(names_core_func, core_at_name_func, staple_at_name_func):
    #Changes the name of the original atoms in the core for AU, and the staples for ST. This names are required for the rest of LAMPIT
    for i in range(len(names_core_func)):
        if(names_core_func[i]==core_at_name_func):
            names_core_func[i]='AU'
        elif (names_core_func[i]==staple_at_name_func):
            names_core_func[i]='ST'
    return names_core_func

def get_anchor(xyz_core_func, names_core_func, name_anchor_tmp):
    #Returns xyz coordinates of the anchors
    return xyz_core_func[names_core_func==name_anchor_tmp,:]

def get_stone(xyz_stones_tmp, xyz_lig_tmp, stone_ndx, stones_ndx_func):
    #Returns xyz coordinates of one stone given the ones of the respective anchor and the stone number
    mag_C=np.linalg.norm(xyz_stones_tmp)
    scaling=(mag_C+np.linalg.norm(xyz_lig_tmp[stones_ndx_func[stone_ndx],:]-xyz_lig_tmp[stones_ndx_func[0],:]))/mag_C
    return scaling*xyz_stones_tmp

def get_stones_mol2(xyz_lig_func, ndx_stones_func):
    xyz_stones_mol2_func=np.zeros((N_stones,3))
    for i in range(len(ndx_stones_func)):
        xyz_stones_mol2_func[i,:]=xyz_lig_func[ndx_stones_func[i],:]
    return xyz_stones_mol2_func

def place_stones(xyz_core_func, names_core_func, xyz_lig_func, names_lig_func, name_anchor_func, stones_ndx_func):
    #Returns the xyz coordinates and names of the ligands in the coating (starting from the anchor)
    xyz_anch=get_anchor(xyz_core_func, names_core_func, name_anchor_func)
    N_stones_func=N_stones*N_S
    xyz_stones_func=np.zeros((N_stones_func,3))
    names_stones_func=np.zeros(N_stones).astype('str')
    for i in range(N_stones):
        names_stones_func[i]=names_lig_func[stones_ndx_func[i]]
    for i in range(N_S):
        xyz_stones_func[i*N_stones,:]=xyz_anch[i,:]
        for j in range(1,N_stones):
            xyz_stones_func[i*N_stones+j,:]=get_stone(xyz_stones_func[i*N_stones,:], xyz_lig_func, j, stones_ndx_func)
    return xyz_stones_func, names_stones_func

def coat_NP(xyz_core_func, names_core_func, xyz_stones_func, xyz_lig_func, names_lig_func, names_stones_func, ndx_stones):
    #Merges xyz coordinates and names of the core and the ligands into one coated NP
    keep_rows=[]
    for i in range(len(names_core_func)):
        if names_core_func[i]=='AU' or names_core_func[i]=='ST':
            keep_rows.append(i)
    xyz_naked_func=xyz_core_func[keep_rows,:]
    names_naked_func=names_core_func[keep_rows]

    xyz_coated_func=xyz_naked_func
    names_coated_func=names_naked_func

    xyz_stones_mol2=get_stones_mol2(xyz_lig_func, ndx_stones)
    xyz_lig_func_conv=np.insert(xyz_lig_func, 3, 1, axis=1).T
    for i in range(N_S):
        xyz_stones_now=xyz_stones_func[i*N_stones:(i+1)*N_stones,:]
        trans_matrix=affine_matrix_from_points(xyz_stones_mol2.T, xyz_stones_now.T, shear=False, scale=False, usesvd=True)
        trans_lig=np.dot(trans_matrix, xyz_lig_func_conv).T[:,:3]

        xyz_coated_func=np.append(xyz_coated_func, trans_lig, axis=0)
        names_coated_func=np.append(names_coated_func, names_lig_func, axis=0)
    return xyz_coated_func, names_coated_func

def print_NP_pdb(xyz_coated_func, names_coated_func, names_stones_func, out_fname):
    #Writes the pdb of the core and the placed stones
    output=open(out_fname, "w")
    res=1
    at=0
    for i in range(len(names_coated_func)):
        at+=1
        at_name_act=names_coated_func[i]
        if at_name_act=='AU' or at_name_act=='ST':
            write_pdb_block(at_name_act, at_name_act, xyz_coated_func[i,:], res, at, out_fname)
            res+=1
        else:
            if(at_name_act==names_stones_func[0]):
                res+=1
            write_pdb_block(at_name_act, res_name_opt, xyz_coated_func[i,:], res, at, out_fname)
    output.close()

def write_pdb_block(atname_func, res_name_func, xyz_func, resnum, atnum, out_filename):
    #Writes the information of one atom in a pdb file
    xyz_func=np.round(xyz_func, decimals=4)
    coords=open(out_filename, 'a')
    coords.write('ATOM'.ljust(6))
    coords.write(str(atnum).rjust(5))
    coords.write('  '+str(atname_func).ljust(3))
    coords.write(' '+str(res_name_func).ljust(3))
    coords.write('  '+str(resnum).rjust(4))
    coords.write('    '+ str(xyz_func[0]).rjust(8))
    coords.write(str(xyz_func[1]).rjust(8))
    coords.write(str(xyz_func[2]).rjust(8)+"\n")
    coords.close()

def print_xyz(coordenadas, nombres, fnombre):
    #Prints a random .xyz file with coordinates and names given
    output=open(fnombre, "w")
    output.write(str(len(nombres)) + "\n \n")
    for i in range(len(nombres)):
        output.write(str(nombres[i]) + " " + str(coordenadas[i,0]) + " " + str(coordenadas[i,1]) + " " + str(coordenadas[i,2]) + "\n")
    output.close()

names_core=change_names(names_core, core_at_name_opt, staple_at_name_opt)
xyz_stones, names_stones=place_stones(xyz_core, names_core, xyz_lig, names_lig, name_anchor_opt, stones_ndx_opt)
xyz_coated_NP, names_coated_NP=coat_NP(xyz_core, names_core, xyz_stones, xyz_lig, names_lig, names_stones, stones_ndx_opt)
#print_xyz(xyz_coated_NP, names_coated_NP, outname_opt)
print_NP_pdb(xyz_coated_NP, names_coated_NP, names_stones, outname_opt)
