import numpy as np
import argparse
import os
import sys
from scipy.spatial import distance
import subunits
import collections

#Imports input file with options
inp=np.genfromtxt(sys.argv[1], dtype="str")

NP_gro_opt = "test/test.gro"
NP_top_opt = "test/test.top"
res_name1_opt = "LF1"
res_name2_opt = "LF2"
ligand1_opt = "test/LF1.mol2"
ligand2_opt = "test/LF2.mol2"
workdir_opt = "test"

def load_options(input_file):
    #Reads the input file with the options (case sensitive) and globally change the values of the variables storing them
    for i in range(len(input_file)):
        if "grofile" in input_file[i]:
            global grofile_opt
            grofile_opt = input_file[i][1]
        if "topfile" in input_file[i]:
            global topfile_opt
            topfile_opt = input_file[i][1]
        if "resname1" in input_file[i]:
            global res_name1_opt
            res_name1_opt = input_file[i][1]
        if "resname2" in input_file[i]:
            global res_name2_opt
            res_name2_opt = input_file[i][1]
        if "ligand1" in input_file[i]:
            global ligand1_opt
            ligand1_opt = input_file[i][1]
        if "ligand2" in input_file[i]:
            global ligand2_opt
            ligand2_opt = input_file[i][1]
        if "workdir" in input_file[i]:
            global workdir_opt
            workdir_opt = input_file[i][1]

load_options(inp)

def angle(a, b, c):
    ba = a - b
    bc = c - b
    cosine_angle = np.dot(ba, bc) / (np.linalg.norm(ba) * np.linalg.norm(bc))
    angle = np.arccos(cosine_angle)
    return angle

def get_lig_info(lig_mol2):
    mol2_file = np.genfromtxt(lig_mol2, delimiter = "\n", dtype='str')
    found_ATOM = 0
    found_CONNECT = 0
    names = []
    types = []
    xyz = []
    for line in mol2_file:
        if "@<TRIPOS>ATOM" in line:
            found_ATOM = 1
        elif "@<TRIPOS>BOND" in line:
            found_ATOM = 0
        elif found_ATOM == 1:
            line = str(line).split()
            xyz.append(line[2:5])
            names.append(line[1])
            types.append(line[5])
        elif "@<TRIPOS>RESIDUECONNECT" in line:
            found_CONNECT = 1
        elif found_CONNECT == 1:
            xyz = np.array(xyz).astype('float')
            ndx_C = int(str(line).split()[0])-1
            dists = distance.cdist([xyz[ndx_C]], xyz)
            ndx_Hs = np.argsort(dists[0])[1:3]
            break
    names = np.array(names)
    types = np.array(types)
    return names[ndx_C], types[ndx_C], names[ndx_Hs], types[ndx_Hs]

def load_gro():
    gro_file = np.genfromtxt(NP_gro_opt, dtype='str', skip_header=2, skip_footer=1)
    xyz = []
    names = []
    resids = []
    for line in gro_file:
        xyz.append(line[4:7])
        names.append(line[2])
        resids.append(line[1])
    return np.array(xyz).astype('float'), np.array(names)

def load_top():
    top_file = np.genfromtxt(NP_top_opt, dtype='str', delimiter = "\n")
    types = []
    for i in range(len(top_file)):
        if "[ atoms ]" in top_file[i]:
            ini = i + 2
        elif "[ bonds ]" in top_file[i]:
            fin = i
    for i in range(ini, fin):
        types.append(top_file[i].split()[1])
    return types

def get_gro_ndx(names_array, search_object):
    ndx = []
    for i in range(len(names_array)):
        if search_object == str(names_array[i]):
            ndx.append(i)
    return np.array(ndx)

def make_blocks(ind_AU, ind_ST, ind_C, ind_H):
    D_ST_AU = distance.cdist(xyz_sys[ind_ST], xyz_sys[ind_AU])
    D_ST_C = distance.cdist(xyz_sys[ind_ST], xyz_sys[ind_C])
    blocks = []
    for i in range(len(ndx_ST)):
        near_Au = D_ST_AU[i].argsort()[0:2]
        near_C = D_ST_C[i].argsort()[0]
        D_C_H = distance.cdist(np.array([xyz_sys[ind_C[near_C]]]), xyz_sys[ind_H])
        near_Hs = D_C_H[0].argsort()[0:2]

        blocks.append(subunits.Block(ndx_S=ind_ST[i], ndx_Au=ind_AU[near_Au], ndx_C=ind_C[near_C], ndx_H=ind_H[near_Hs]))
    return blocks

def make_staples(blocks_list):
    N_blocks = len(blocks_list)
    staples_ndx = []
    for i in range(N_blocks):
        ndx_tmp = []
        for j in range(N_blocks):
            if not set(blocks_list[i].Au).isdisjoint(blocks_list[j].Au):
                ndx_tmp.append(j)
        staples_ndx.append((ndx_tmp))
    staples_ndx = [list(x) for x in set(tuple(x) for x in staples_ndx)] #Takes the unique lists within the list

    staples = []
    for i in range(len(staples_ndx)):
        blocks_S = np.array([])
        blocks_Au = np.array([])
        blocks_C = np.array([])
        blocks_H = np.array([])
        for j in range(len(staples_ndx[i])):
            block_act = blocks_list[staples_ndx[i][j]]
            blocks_S = np.append(blocks_S, block_act.S)
            blocks_Au = np.append(blocks_Au, block_act.Au)
            blocks_C = np.append(blocks_C, block_act.C)
            blocks_H = np.append(blocks_H, block_act.H)
        blocks_Au_l = [item for item, count in collections.Counter(blocks_Au).items() if count > 1]
        blocks_S = np.unique(blocks_S)
        blocks_Au = np.unique(blocks_Au)
        blocks_C = np.unique(blocks_C)
        blocks_H = np.unique(blocks_H)

        staples.append(subunits.Staple(tipo='UNK', ndx_S=blocks_S, ndx_Au=blocks_Au, ndx_Au_l=blocks_Au_l, ndx_C=blocks_C, ndx_H=blocks_H))
    return staples

def classify_staples(staples_list):
    for i in range(len(staples_list)):
        staple_act = staples_list[i]
        N_S = len(staple_act.S)
        N_Au = len(staple_act.Au)
        if N_S == 1 and N_Au == 2:
            staple_act.change_tipo('STP')
        elif N_S == 2 and N_Au == 3:
            staple_act.change_tipo('STR')
        elif N_S == 3 and N_Au == 4:
            D_S_Aul = distance.cdist(xyz_sys[staple_act.S], xyz_sys[staple_act.Au_l])
            for j in range (len(staple_act.S)):
                near_Au = staple_act.S[D_Aul_S[j].argsort()[0:2]]
                if np.all(np.in1d(near_Au, staple_act.Au_l)):
                    angle = angle(xyz_sys[near_Au[0]], xyz_sys[staple_act.S[j]], xyz_sys[near_Au[1]]))
                    if angle <= 109.0 and angle >= 91.0:
                        staple_act.change_tipo('STC')
                    elif angle <= 128.2 and angle >= 110.2:
                        staple_act.change_tipo('STV')
                    else:
                        print("Unrecognized staple")
        else:
            print("Unrecognized staple")
    return staples_list

def staple_to_residues(staples_list):
    residues = []
    for i in range(len(staples_list)):
        sta = staples_list[i]
        restype = sta.tipo
        ndx = np.concatenate((sta.S, sta.Au, sta.C, sta.H)).astype('int')
        residues.append(subunits.Residue(restype=restype, ndx=ndx))
    return residues

def print_pdb(residues_list, out_fname):
    #Writes a pdb file with only the staples and each one with their resID and restype
    at=0
    resnum = 0
    for i in range(len(residues_list)):
        res_act = residues_list[i]
        ndx = res_act.ndx
        resnum += 1
        for j in range(len(ndx)):
            at+=1
            write_pdb_block(names_sys[ndx[j]], res_act.restype, xyz_sys[ndx[j]], resnum, at, out_fname)
    output=open(out_fname, 'a')
    output.write('END')
    output.close()

def write_pdb_block(atname_func, res_name_func, xyz_func, resnum, atnum, out_filename):
    #Writes one line of a generic pdb file
    xyz_func=np.round(xyz_func*10., decimals=4)
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

def write_bonds(staples_list):
    bonds = open(workdir_opt + '/bonds.top', 'w')
    func_type = str(1)
    for i in range(len(staples_list)):
        s = staples_list[i]
        D_S_Au = distance.cdist(xyz_sys[s.S], xyz_sys[s.Au])
        D_S_C = distance.cdist(xyz_sys[s.S], xyz_sys[s.C])
        for j in range(len(s.S)):
            near_Au = s.Au[D_S_Au[j].argsort()[0:2]]
            cons = 62730
            for k in range(len(near_Au)):
                if near_Au[k] in s.Au_l:
                    zero = 0.233
                else:
                    zero = 0.241
                bonds.write(str(s.S[j]+1).rjust(6)+str(near_Au[k]+1).rjust(7)+str(func_type).rjust(4)+"{:.4e}".format(zero).rjust(14)+"{:.4e}".format(cons).rjust(14)+" ; "+names_sys[s.S[j]]+" - "+names_sys[near_Au[k]]+" m\n")

            near_C = s.C[D_S_C[j].argsort()[0]]
            if types_sys[near_C] == 'CT':
                cons = 99113.0
                zero = 0.184
                bonds.write(str(s.S[j]+1).rjust(6)+str(near_C+1).rjust(7)+str(func_type).rjust(4)+"{:.4e}".format(zero).rjust(14)+"{:.4e}".format(cons).rjust(14)+" ; "+names_sys[s.S[j]]+" - "+names_sys[near_C]+" m\n")

def write_angles(staples_list):
    angles = open(workdir_opt + '/angles.top', 'w')
    func_type = str(1)
    for i in range(len(staples_list)):
        s = staples_list[i]
        D_S_Au = distance.cdist(xyz_sys[s.S], xyz_sys[s.Au])
        D_S_C = distance.cdist(xyz_sys[s.S], xyz_sys[s.C])
        for j in range(len(s.S)):
            #Au-S-Au
            near_Au = s.Au[D_S_Au[j].argsort()[0:2]]
            cons = 460.24
            if ((near_Au[0] in s.Au_l) and (near_Au[1] not in s.Au_l) or (near_Au[0] not in s.Au_l) and (near_Au[1] in s.Au_l)):
                zero = 91.3
            elif ((near_Au[0] in s.Au_l) and (near_Au[1] in s.Au_l)):
                if s.tipo == "STV":
                    cons = 1460.24
                    zero = 119.2
                elif s.tipo == "STC":
                    zero = 100.0
            else:
                print("There is an unsupported Au-S-Au bond")
            angles.write(str(near_Au[0]+1).rjust(6)+str(s.S[j]+1).rjust(7)+str(near_Au[1]+1).rjust(7)+str(func_type).rjust(7)+"{:.4e}".format(zero).rjust(14)+"{:.4e}".format(cons).rjust(14)+" ; "+names_sys[near_Au[0]]+" - "+names_sys[s.S[j]]+" - "+names_sys[near_Au[1]]+" m\n")

            near_C = s.C[D_S_C[j].argsort()[0]]
            cons = 146.37
            for k in range(len(near_Au)):
                if near_Au[k] in s.Au_l:
                    zero = 106.8
                else:
                    zero = 111.6
                angles.write(str(near_Au[k]+1).rjust(6)+str(s.S[j]+1).rjust(7)+str(near_C+1).rjust(7)+str(func_type).rjust(7)+"{:.4e}".format(zero).rjust(14)+"{:.4e}".format(cons).rjust(14)+" ; "+names_sys[near_Au[k]]+" - "+names_sys[s.S[j]]+" - "+names_sys[near_C]+" m\n")

            D_C_H = distance.cdist([xyz_sys[near_C]], xyz_sys[s.H])
            near_H = s.H[D_C_H[0].argsort()[0:2]]
            cons = 418.40
            zero = 107.0
            for k in range(len(near_H)):
                if types_sys[near_H[k]] == "HC":
                    angles.write(str(s.S[j]+1).rjust(6)+str(near_C+1).rjust(7)+str(near_H[k]+1).rjust(7)+str(func_type).rjust(7)+"{:.4e}".format(zero).rjust(14)+"{:.4e}".format(cons).rjust(14)+" ; "+names_sys[s.S[j]]+" - "+names_sys[near_C]+" - "+names_sys[near_H[k]]+" m\n")
                else:
                    print("There are no parameters for angles involving this kind of hydrogen atoms")

        D_Aul_S = distance.cdist(xyz_sys[s.Au_l], xyz_sys[s.S])
        cons = 460.240
        zero = 172.4
        for j in range(len(s.Au_l)):
            near_S = s.S[D_Aul_S[j].argsort()[0:2]]
            angles.write(str(near_S[0]+1).rjust(6)+str(s.Au_l[j]+1).rjust(7)+str(near_S[1]+1).rjust(7)+str(func_type).rjust(7)+"{:.4e}".format(zero).rjust(14)+"{:.4e}".format(cons).rjust(14)+" ; "+names_sys[near_S[0]]+" - "+names_sys[s.Au_l[j]]+" - "+names_sys[near_S[1]]+" m\n")

def write_topology():

    top_file = open(NP_top_opt, "r")
    cont_top_file = top_file.readlines()
    top_file.close()

    final_top=open(NP_top_opt, "w")
    final_top.close()
    final_top=open(NP_top_opt, "a")

    for i in range(len(cont_top_file)):
        final_top.writelines(cont_top_file[i])
        if ";   ai     aj funct   r             k" in cont_top_file[i]:
            bonds_file=open(workdir_opt + "/bonds.top","r")
            bonds_contents=bonds_file.readlines()
            bonds_file.close()
            for j in range(len(bonds_contents)):
                final_top.writelines(bonds_contents[j])
        if ";   ai     aj     ak    funct   theta         cth" in cont_top_file[i]:
            angles_file=open(workdir_opt + "/angles.top","r")
            angles_contents=angles_file.readlines()
            angles_file.close()
            for j in range(len(angles_contents)):
                final_top.writelines(angles_contents[j])


name_C1 , type_C1, name_H1, type_H1 = get_lig_info(ligand1_opt)
name_C2 , type_C2, name_H2, type_H2 = get_lig_info(ligand2_opt)

xyz_sys, names_sys = load_gro()
types_sys = load_top()

ndx_AU = get_gro_ndx(names_sys, 'AU')
ndx_ST = get_gro_ndx(names_sys, 'ST')
ndx_C1 = get_gro_ndx(names_sys, name_C1)
ndx_C2 = get_gro_ndx(names_sys, name_C2)
ndx_H1a = get_gro_ndx(names_sys, name_H1[0])
ndx_H1b = get_gro_ndx(names_sys, name_H1[1])
ndx_H2a = get_gro_ndx(names_sys, name_H2[0])
ndx_H2b = get_gro_ndx(names_sys, name_H2[1])

ndx_C = np.unique(np.append(ndx_C1, ndx_C2))
ndx_H = np.unique(np.append(ndx_H1a, [ndx_H1b, ndx_H2a, ndx_H2b]))

blocks = make_blocks(ndx_AU, ndx_ST, ndx_C, ndx_H)
staples = make_staples(blocks)
staples = classify_staples(staples)
residues = staple_to_residues(staples)
print_pdb(residues, workdir_opt + '/staples.pdb')

write_bonds(staples)
write_angles(staples)
write_topology()
