import numpy as np
import matplotlib.pyplot as plt
import argparse
import os
import sys
from  transformations import *
import random
from mpl_toolkits.mplot3d import Axes3D

#Declaration of the flags to run the script

inp=np.genfromtxt(sys.argv[1], dtype="str")
outname_opt = "Output.pdb"
corename_opt = "Core.pdb"
core_at_name_opt = "Au"
staple_at_name_opt = "S"
name_anchor_opt = "C"
res_name1_opt = "UNK"
res_name2_opt = "UNK"
ligname1_opt = "Ligand1.mol2"
ligname2_opt = "Ligand2.mol2"
frac_lig1_opt = 0.5
morph_opt = "random"
rseed_opt = 666

def load_options(input_file):
    for i in range(len(input_file)):
        if "output" in input_file[i]:
            global outname_opt
            outname_opt = input_file[i][1]
        elif "core" in input_file[i]:
            global corename_opt
            corename_opt = input_file[i][1]
        elif "corename" in input_file[i]:
            global core_at_name_opt
            core_at_name_opt = input_file[i][1]
        elif "staplename" in input_file[i]:
            global staple_at_name_opt
            staple_at_name_opt = input_file[i][1]
        elif "coreanchor" in input_file[i]:
            global name_anchor_opt
            name_anchor_opt = input_file[i][1]
        elif "ligname1" in input_file[i]:
            global res_name1_opt
            res_name1_opt = input_file[i][1]
        elif "ligname2" in input_file[i]:
            global res_name2_opt
            res_name2_opt = input_file[i][1]
        elif "ligand1" in input_file[i]:
            global ligname1_opt
            ligname1_opt = input_file[i][1]
        elif "ligand2" in input_file[i]:
            global ligname2_opt
            ligname2_opt = input_file[i][1]
        elif "morphology" in input_file[i]:
            global morph_opt
            morph_opt = input_file[i][1]
        elif "frac_lig1" in input_file[i]:
            global frac_lig1_opt
            frac_lig1_opt = float(input_file[i][1])
        elif "rseed" in input_file[i]:
            global rseed_opt
            rseed_opt = int(input_file[i][1])

load_options(inp)

def check_options(core_fname, ligand_fname1, ligand_fname2):
    #Checks that the required options exist and are consistent
    if not (os.path.isfile(core_fname)):
        sys.exit("Core file could not be found")
    if not (os.path.isfile(ligand_fname1)):
        sys.exit("Ligand file 1 could not be found")
    if not (os.path.isfile(ligand_fname2)):
        sys.exit("Ligand file 2 could not be found")

def change_names(names_core_func):
    #Changes the name of the original atoms in the core for AU, and the staples for ST. This names are required for the rest of LAMPIT
    for i in range(len(names_core_func)):
        if(names_core_func[i]==core_at_name_opt):
            names_core_func[i]='AU'
        elif (names_core_func[i]==staple_at_name_opt):
            names_core_func[i]='ST'
    return names_core_func

def init_lig_mol2(ligand_fname):
    #Imports ligand mol2 file. Returns xyz coordinates, names, and index corresponding to the anchor
    lig_file=np.genfromtxt(ligand_fname, delimiter='\n', dtype='str')
    N_lig_file=len(lig_file)
    found_ATOM=0
    names_lig_func=[]
    xyz_lig_func=[]
    for i in range(N_lig_file):
        if "@<TRIPOS>BOND" in lig_file[i]:
            fin = i
            found_ATOM = 0
        elif found_ATOM == 1:
            at_file=lig_file[i].split()
            names_lig_func.append(at_file[1])
            xyz_lig_func.append(at_file[2:5])
        elif "@<TRIPOS>ATOM" in lig_file[i]:
            ini = i+1
            found_ATOM = 1
        elif "@<TRIPOS>RESIDUECONNECT" in lig_file[i]:
            anchor_ndx_func=int(lig_file[i+1].split()[0])-1

    xyz_lig_func, names_lig_func=np.array(xyz_lig_func, dtype='float'), np.array(names_lig_func)
    anchor_pos=np.copy(xyz_lig_func)[anchor_ndx_func,:]

    #Moves the ligand so that the anchor is in (0,0,0)
    for i in range(len(xyz_lig_func[:,0])):
        xyz_lig_func[i,:] = xyz_lig_func[i,:] - anchor_pos
    return xyz_lig_func, names_lig_func, anchor_ndx_func

def init_core_pdb(core_fname):
    #Imports core pdb file. Centers the core in (0,0,0) and returns xyz coordinates and names
    core_file=np.genfromtxt(core_fname, delimiter='\n', dtype=str)
    first=1
    N_core_func=0
    for i in range(len(core_file)):
        if "ATOM" in core_file[i] or "HETATM" in core_file[i]:
            N_core_func+=1
            if first==1:
                ini=i
                first=0
    core_file=core_file[ini:  ini+N_core_func]
    names_core_func=[]
    xyz_core_func=np.zeros((N_core_func,3))
    for i in range(N_core_func):
        at_file=core_file[i].split()
        names_core_func.append(at_file[2])
        xyz_core_func[i,:]=at_file[5:8]
    xyz_core_func=np.array(xyz_core_func)
    names_core_func=np.array(names_core_func)

    names_core_func=change_names(names_core_func)

    COM=np.average(xyz_core_func[names_core_func=='AU',:], axis=0)

    for i in range(len(xyz_core_func[:,0])):
        xyz_core_func[i,:]=xyz_core_func[i,:]-COM
    return xyz_core_func, names_core_func

def get_ligand_pill(xyz_lig_func, anchor_ndx_func):
    #Runs a PCA and takes the first eigenvector as the best fitting line.
    pca = np.linalg.eig(np.cov(xyz_lig_func.T))
    pca1 = pca[1][0]
    var1 = pca[0][0]/np.sum(pca[0])*100
    print("The PCA1 explains: {:.1f}% of the points variance".format(var1))

    #Randomly takes 2 other atoms in the ligand and project their positions in PCA1
    random.seed(666)
    rango = list(range(len(xyz_lig_func[:,0])))
    rango.remove(anchor_ndx_func)
    pillars_ndx = random.sample(rango, 2)
    pillars_func = np.array([0.0, 0.0, 0.0]) #This corresponds to the first stone (i.e. the anchor) which will always be in (0,0,0)
    for i in pillars_ndx:
        pillars_func = np.vstack((pillars_func, np.dot(xyz_lig_func[i], pca1) * pca1))
    return pillars_func

check_options(corename_opt, ligname1_opt, ligname2_opt)

xyz_lig1, names_lig1, anchor_ndx1 = init_lig_mol2(ligname1_opt)
xyz_lig2, names_lig2, anchor_ndx2 = init_lig_mol2(ligname2_opt)

xyz_core, names_core = init_core_pdb(corename_opt)

xyz_pillars1 = get_ligand_pill(xyz_lig1, anchor_ndx1)
xyz_pillars2 = get_ligand_pill(xyz_lig2, anchor_ndx2)

#Define number of atoms in a ligand, number of staples, and number of stones specified
N_at_lig1 = len(xyz_lig1[:,0])
N_at_lig2 = len(xyz_lig2[:,0])
N_S = len(names_core[names_core=='ST'])

def assign_morph(xyz_core_func, names_core_func):
    xyz_anchors_func = xyz_core_func[names_core_func==name_anchor_opt,:]
    N_anchors = len(xyz_anchors_func)
    for_lig1 = round(N_anchors*frac_lig1_opt)
    indexes = list(range(N_anchors))
    if morph_opt == "random":
        random.seed(rseed_opt)
        random.shuffle(indexes)
        lig1_ndx = indexes[:for_lig1]
        lig2_ndx = indexes[for_lig1:]
        xyz_anchors1_func=xyz_anchors_func[lig1_ndx,:]
        xyz_anchors2_func=xyz_anchors_func[lig2_ndx,:]
    else:
        print("Morphology not supported")
    return xyz_anchors1_func, xyz_anchors2_func

def get_stones(xyz_anchorsi_func, xyz_pillarsi_func):
    #Return a 3D array with the xyz coordinates for all the stones of all the anchors
    n_stones_lig = len(xyz_pillarsi_func)
    n_anchors = len(xyz_anchorsi_func)
    xyz_stones = np.zeros((n_anchors, n_stones_lig, 3))

    #Takes the COM-C vectors and scale them to match the distance between staples in the ligand's file
    for i in range(n_anchors):
        mag_C = np.linalg.norm(xyz_anchorsi_func[i,:])
        for j in range(n_stones_lig):
            scaling = (mag_C + np.linalg.norm(xyz_pillarsi_func[j,:]))/mag_C
            xyz_stones[i,j,:]=xyz_anchorsi_func[i,:]*scaling
    return xyz_stones

def coat_NP(xyz_core_func, names_core_func, xyz_lig1_func, names_lig1_func, xyz_pillars1_func, xyz_stones1_func, xyz_lig2_func, names_lig2_func, xyz_pillars2_func, xyz_stones2_func):
    #Merges xyz coordinates and names of the core and the ligands into one coated NP
    keep_rows=[]
    for i in range(len(names_core_func)):
        if names_core_func[i]=='AU' or names_core_func[i]=='ST':
            keep_rows.append(i)

    xyz_coated_func=xyz_core_func[keep_rows,:]
    names_coated_func=names_core_func[keep_rows]

    xyz_lig1_func_conv=np.insert(xyz_lig1_func, 3, 1, axis=1).T
    for i in range(len(xyz_stones1_func[:,0,0])):
        xyz_stones_now = xyz_stones1_func[i,:,:]
        trans_matrix=affine_matrix_from_points(xyz_pillars1_func.T, xyz_stones_now.T, shear=False, scale=False, usesvd=True)
        trans_lig=np.dot(trans_matrix, xyz_lig1_func_conv).T[:,:3]

        xyz_coated_func=np.append(xyz_coated_func, trans_lig, axis=0)
        names_coated_func=np.append(names_coated_func, names_lig1_func, axis=0)

    xyz_lig2_func_conv=np.insert(xyz_lig2_func, 3, 1, axis=1).T
    for i in range(len(xyz_stones2_func[:,0,0])):
        xyz_stones_now = xyz_stones2_func[i,:,:]
        trans_matrix=affine_matrix_from_points(xyz_pillars2_func.T, xyz_stones_now.T, shear=False, scale=False, usesvd=True)
        trans_lig=np.dot(trans_matrix, xyz_lig2_func_conv).T[:,:3]

        xyz_coated_func=np.append(xyz_coated_func, trans_lig, axis=0)
        names_coated_func=np.append(names_coated_func, names_lig2_func, axis=0)
    return xyz_coated_func, names_coated_func

def print_NP_pdb(xyz_coated_func, names_coated_func, xyz_anchors1_func, xyz_anchors2_func, out_fname):
    #Writes the pdb of the core and the placed stones
    N_lig1 = len(xyz_anchors1)
    N_tot_lig1 = N_lig1 * N_at_lig1
    N_lig2 = len(xyz_anchors2)
    N_tot_lig2 = N_lig2 * N_at_lig2
    N_core = len(xyz_coated_func) - N_lig1*N_at_lig1 - N_lig2*N_at_lig2

    res=0
    at=0
    output=open(out_fname, "w")
    for i in range(N_core):
        at+=1
        res+=1
        at_name_act=names_coated_func[i]
        write_pdb_block(at_name_act, at_name_act, xyz_coated_func[i,:], res, at, out_fname)

    lig_atoms=0
    for i in range(N_tot_lig1):
        at+=1
        at_name_act=names_coated_func[i+N_core]
        if(lig_atoms%N_at_lig1==0):
            res+=1
        lig_atoms+=1
        write_pdb_block(at_name_act, res_name1_opt, xyz_coated_func[i+N_core,:], res, at, out_fname)

    lig_atoms=0
    for i in range(N_tot_lig2):
        at+=1
        at_name_act=names_coated_func[i+N_core+N_tot_lig1]
        if(lig_atoms%N_at_lig2==0):
            res+=1
        lig_atoms+=1
        write_pdb_block(at_name_act, res_name2_opt, xyz_coated_func[i+N_core+N_tot_lig1,:], res, at, out_fname)

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
    coords.write('    '+ str(round(xyz_func[0],3)).rjust(8))
    coords.write(str(round(xyz_func[1],3)).rjust(8))
    coords.write(str(round(xyz_func[2],3)).rjust(8)+"\n")
    coords.close()

def print_xyz(coordenadas, nombres, fnombre):
    #Prints a random .xyz file with coordinates and names given
    output=open(fnombre, "w")
    output.write(str(len(nombres)) + "\n \n")
    for i in range(len(nombres)):
        output.write(str(nombres[i]) + " " + str(coordenadas[i,0]) + " " + str(coordenadas[i,1]) + " " + str(coordenadas[i,2]) + "\n")
    output.close()

xyz_anchors1, xyz_anchors2 = assign_morph(xyz_core, names_core)

xyz_stones1 = get_stones(xyz_anchors1, xyz_pillars1)
xyz_stones2 = get_stones(xyz_anchors2, xyz_pillars2)

xyz_coated_NP, names_coated_NP = coat_NP(xyz_core, names_core, xyz_lig1, names_lig1, xyz_pillars1, xyz_stones1, xyz_lig2, names_lig2, xyz_pillars2, xyz_stones2)

print_NP_pdb(xyz_coated_NP, names_coated_NP, xyz_anchors1, xyz_anchors2, outname_opt)
#print_xyz(xyz_coated_NP, names_coated_NP, outname_opt)
