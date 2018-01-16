#imports necessary packages and classes
import numpy as np
from optparse import OptionParser
from scipy.spatial import distance
import atom
import matplotlib.pyplot as plt
import os

#Declaration of the flags to run the script
parser=OptionParser()
parser.add_option("-x", "--gro", action="store", type='string', dest="GroFile", default='NP.gro', help="Name of the .gro file with the staples")
parser.add_option("-p", "--top", action="store", type='string', dest="TopFile", default='NP.top', help="Name of the .top file associated to the file in -x")
parser.add_option("-a", "--anchor", action="store", type='string', dest="AnchorName", default='C1', help="Name of the ligands' atom linked to the staple")
parser.add_option("-f", "--folder", action="store", type='string', dest="WorkDir", default='NP', help="Path of the working directory to save intermediate files")
parser.add_option("-y", "--hydrogen", action="append", dest="ListHAnchor", default=[], help="List of names of hydrogen atoms involved in a 3-atom angle with an ST")
(options, args)= parser.parse_args()
grofile_opt=options.GroFile
topfile_opt=options.TopFile
anchorname_opt=options.AnchorName
workdir_opt=options.WorkDir
list_Hs_opt=options.ListHAnchor[0].split(',')

def init_element(grofile_func, element):
    #Returns the xyz coordinates of the atoms with name 'element'
    elem_file_func=np.genfromtxt(grofile_func, delimiter='\n', dtype='str', skip_header=2)
    N_elem_func=0
    ndx_element=[]
    for i in range(len(elem_file_func)):
        line=elem_file_func[i].split()
        if element in line:
            N_elem_func+=1
            ndx_element.append(i)
    elem_xyz_func=np.zeros((N_elem_func,3))
    for i in range(len(ndx_element)):
        elem_actual=elem_file_func[ndx_element[i]].split()
        elem_xyz_func[i,:]=elem_actual[4:7]

    return np.array(elem_xyz_func)*10

def init_topology(topfile_func):
    #Returns two lists with all the names and types of the atoms in the system in the same order as in the gro fiel
    names_file_func=np.genfromtxt(topfile_func, delimiter='\n', dtype='str')
    for i in range(len(names_file_func)):
        if "[ atoms ]" in names_file_func[i]:
            ini=i+2
        elif "[ bonds ]" in names_file_func[i]:
            fin=i
    names_all_func=[]
    types_all_func=[]
    for i in range(ini, fin):
        names_all_func.append(names_file_func[i].split()[4])
        types_all_func.append(names_file_func[i].split()[1])
    return np.array(names_all_func), np.array(types_all_func)

def calculate_angles(pa, pb, pc):
    #Given 3 arrays Nx3, it returns an Nx1 array with the angle made by pa[i], pb[i], and pc[i], centered in pa[i]
    N_angles=len(pa[:,0])
    ab=np.subtract(pb, pa)
    ac=np.subtract(pc,pa)
    cos_ang=np.zeros(N_angles)
    for i in range(N_angles):
        cos_ang[i]=np.dot(ab[i,:], ac[i,:])/(np.linalg.norm(ab[i,:])*np.linalg.norm(ac[i,:]))
    return np.arccos(cos_ang)*180/np.pi

def calculate_distance(pa, pb, pc):
    #Given 3 arrays Nx3, an Nx2 array is returned with the distance between pa[i]-pb[i] (col 1) and pa[i]-pc[i] (col 2)
    N_dist=len(pa[:,0])
    dist_func=np.zeros((N_dist,2))
    pab=np.subtract(pb,pa)
    pac=np.subtract(pc,pa)
    for i in range(N_dist):
        dist_func[i,0]=np.linalg.norm(pab[i,:])
        dist_func[i,1]=np.linalg.norm(pac[i,:])
    return dist_func

def print_staple_pdb(residues_list, AU_xyz_func, ST_xyz_func, out_fname):
    #Writes a pdb file with only the staples and each one with their resID and restype
    at=0
    for i in range(len(residues_list)):
        res_act=residues_list[i]
        for j in range(len(res_act.au_atoms)):
            at+=1
            write_pdb_block('AU', res_act.restype, AU_xyz_func[res_act.au_atoms[j],:], res_act.resnum, at, out_fname)
        for j in range(len(res_act.s_atoms)):
            at+=1
            write_pdb_block('ST', res_act.restype, ST_xyz_func[res_act.s_atoms[j],:], res_act.resnum, at, out_fname)
    output=open(out_fname, 'a')
    output.write('END')
    output.close()

def write_pdb_block(atname_func, res_name_func, xyz_func, resnum, atnum, out_filename):
    #Writes one line of a generic pdb file
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

def print_new_bonds(res_list_tmp, all_mins_tmp, mins_anchor_tmp, bonds_file_func, names_all_func, types_all_func, anch_name):
    #Prints the bonds with staple atoms
    bonds_func=open(bonds_file_func, 'w')
    function_type=str(1)
    for res_act in res_list_tmp:
        for j in range(len(res_act.s_atoms)):
            s_atom_act=res_act.s_atoms[j]
            real_s_ndx=np.where(names_all_func=='ST')[0][s_atom_act]+1
            for k in range(2):
                au_ndx = all_mins_tmp[s_atom_act,k]
                real_au_ndx=np.where(names_all_func=='AU')[0][au_ndx]+1
                cons_value=62730
                if res_act.au_n_bonds[np.where(res_act.au_atoms==au_ndx)]==2:
                    zero_value=0.233
                elif res_act.au_n_bonds[np.where(res_act.au_atoms==au_ndx)]==1:
                    zero_value=0.241
                else:
                    print('WTF. Gold atoms with unknown valence')
                bonds_func.write(str(real_s_ndx).rjust(6)+str(real_au_ndx).rjust(7)+str(function_type).rjust(4)+"{:.4e}".format(zero_value).rjust(14)+"{:.4e}".format(cons_value).rjust(14)+" ; "+names_all_func[real_s_ndx-1]+" - "+names_all_func[real_au_ndx-1]+" m\n")

            anch_ndx=mins_anchor_tmp[s_atom_act]
            real_anchor_ndx=np.where(names_all_func==anch_name)[0][anch_ndx]+1
            if types_all_func[np.where(names_all_func==anch_name)[0][anch_ndx]]=='CA':
                cons_value=198321.6
                zero_value=0.175
            elif types_all_func[np.where(names_all_func==anch_name)[0][anch_ndx]]=='CT':
                cons_value=99113.0
                zero_value=0.184
            else:
                print('Missing parameters for ST-anchor bond')
            bonds_func.write(str(real_s_ndx).rjust(6)+str(real_anchor_ndx).rjust(7)+str(function_type).rjust(4)+"{:.4e}".format(zero_value).rjust(14)+"{:.4e}".format(cons_value).rjust(14)+" ; "+names_all_func[real_s_ndx-1]+" - "+names_all_func[real_anchor_ndx-1]+" m\n")
    bonds_func.close()

def print_new_angles(res_list_tmp, all_mins_tmp, mins_anchor_tmp, angles_file_func, names_all_func, types_all_func, anch_name, list_H_names_func):
    #Prints the angles with staple atoms present
    angles_func=open(angles_file_func, "w")
    function_type=str(1)
    for res_act in res_list_tmp:
        saus_list=[]
        for j in range(len(res_act.s_atoms)):
            s_atom_act=res_act.s_atoms[j]
            real_s_ndx=np.where(names_all_func=='ST')[0][s_atom_act]+1
            anch_ndx=mins_anchor_tmp[s_atom_act]
            real_anchor_ndx=np.where(names_all_func==anch_name)[0][anch_ndx]+1

            au_ndx1=all_mins_tmp[s_atom_act,0]
            au_ndx2=all_mins_tmp[s_atom_act,1]
            neigh_count_au1=res_act.au_n_bonds[np.where(res_act.au_atoms==au_ndx1)]
            neigh_count_au2=res_act.au_n_bonds[np.where(res_act.au_atoms==au_ndx2)]
            real_au_ndx1=np.where(names_all_func=='AU')[0][au_ndx1]+1
            real_au_ndx2=np.where(names_all_func=='AU')[0][au_ndx2]+1

            #Au-S-AU angles
            if (neigh_count_au1==2) and (neigh_count_au2==2):
                if res_act.restype == 'STC':
                    cons_value=460.240
                    zero_value=100.0
                elif res_act.restype == 'STV':
                    cons_value=1460.240
                    zero_value=119.2
            elif (((neigh_count_au1==1) and (neigh_count_au2==2)) or ((neigh_count_au1==2) and (neigh_count_au2==1))):
                cons_value=460.240
                zero_value=91.3
            elif ((neigh_count_au1==1) and (neigh_count_au2==1)):
                print('There are unparametrized STP residues')
            angles_func.write(str(real_au_ndx1).rjust(6)+str(real_s_ndx).rjust(7)+str(real_au_ndx2).rjust(7)+str(function_type).rjust(7)+"{:.4e}".format(zero_value).rjust(14)+"{:.4e}".format(cons_value).rjust(14)+" ; "+names_all_func[real_au_ndx1-1]+" - "+names_all_func[real_s_ndx-1]+" - "+names_all_func[real_au_ndx2-1]+" m\n")

            #S-AU-S angles
            if (neigh_count_au1==2):
                s_atoms_angle=np.where(all_mins_tmp==au_ndx1)[0]
                real_s_ndx1=np.where(names_all_func=='ST')[0][s_atoms_angle[0]]+1
                real_s_ndx2=np.where(names_all_func=='ST')[0][s_atoms_angle[1]]+1
                saus_sum=real_s_ndx1+real_au_ndx1+real_s_ndx2
                if saus_sum not in saus_list:
                    saus_list.append(saus_sum)
                    cons_value=460.240
                    zero_value=172.4
                    angles_func.write(str(real_s_ndx1).rjust(6)+str(real_au_ndx1).rjust(7)+str(real_s_ndx2).rjust(7)+str(function_type).rjust(7)+"{:.4e}".format(zero_value).rjust(14)+"{:.4e}".format(cons_value).rjust(14)+" ; "+names_all_func[real_s_ndx1-1]+" - "+names_all_func[real_au_ndx1-1]+" - "+names_all_func[real_s_ndx2-1]+" m\n")

                #AU_Lig-S-C
                cons_value=146.370
                zero_value=106.8
                angles_func.write(str(real_au_ndx1).rjust(6)+str(real_s_ndx).rjust(7)+str(real_anchor_ndx).rjust(7)+str(function_type).rjust(7)+"{:.4e}".format(zero_value).rjust(14)+"{:.4e}".format(cons_value).rjust(14)+" ; "+names_all_func[real_au_ndx1-1]+" - "+names_all_func[real_s_ndx-1]+" - "+names_all_func[real_anchor_ndx-1]+" m\n")

            if (neigh_count_au2==2):
                s_atoms_angle=np.where(all_mins_tmp==au_ndx2)[0]
                real_s_ndx1=np.where(names_all_func=='ST')[0][s_atoms_angle[0]]+1
                real_s_ndx2=np.where(names_all_func=='ST')[0][s_atoms_angle[1]]+1
                saus_sum=real_s_ndx1+real_au_ndx2+real_s_ndx2
                if saus_sum not in saus_list:
                    saus_list.append(saus_sum)
                    cons_value=460.240
                    zero_value=172.4
                    angles_func.write(str(real_s_ndx1).rjust(6)+str(real_au_ndx2).rjust(7)+str(real_s_ndx2).rjust(7)+str(function_type).rjust(7)+"{:.4e}".format(zero_value).rjust(14)+"{:.4e}".format(cons_value).rjust(14)+" ; "+names_all_func[real_s_ndx1-1]+" - "+names_all_func[real_au_ndx2-1]+" - "+names_all_func[real_s_ndx2-1]+" m\n")

                #AU_Lig-S-C
                cons_value=146.370
                zero_value=106.8
                angles_func.write(str(real_au_ndx2).rjust(6)+str(real_s_ndx).rjust(7)+str(real_anchor_ndx).rjust(7)+str(function_type).rjust(7)+"{:.4e}".format(zero_value).rjust(14)+"{:.4e}".format(cons_value).rjust(14)+" ; "+names_all_func[real_au_ndx2-1]+" - "+names_all_func[real_s_ndx-1]+" - "+names_all_func[real_anchor_ndx-1]+" m\n")

            #Au_surf-S-C
            if (neigh_count_au1==1):
                cons_value=146.370
                zero_value=111.6
                angles_func.write(str(real_au_ndx1).rjust(6)+str(real_s_ndx).rjust(7)+str(real_anchor_ndx).rjust(7)+str(function_type).rjust(7)+"{:.4e}".format(zero_value).rjust(14)+"{:.4e}".format(cons_value).rjust(14)+" ; "+names_all_func[real_au_ndx1-1]+" - "+names_all_func[real_s_ndx-1]+" - "+names_all_func[real_anchor_ndx-1]+" m\n")
            if (neigh_count_au2==1):
                cons_value=146.370
                zero_value=111.6
                angles_func.write(str(real_au_ndx2).rjust(6)+str(real_s_ndx).rjust(7)+str(real_anchor_ndx).rjust(7)+str(function_type).rjust(7)+"{:.4e}".format(zero_value).rjust(14)+"{:.4e}".format(cons_value).rjust(14)+" ; "+names_all_func[real_au_ndx2-1]+" - "+names_all_func[real_s_ndx-1]+" - "+names_all_func[real_anchor_ndx-1]+" m\n")

            #S-CT-HC
            for an_H in list_H_names_func:
                cons_value=418.40
                zero_value=107.0
                real_H_ndx=np.where(names_all_func==an_H)[0][anch_ndx]+1
                angles_func.write(str(real_s_ndx).rjust(6)+str(real_anchor_ndx).rjust(7)+str(real_H_ndx).rjust(7)+str(function_type).rjust(7)+"{:.4e}".format(zero_value).rjust(14)+"{:.4e}".format(cons_value).rjust(14)+" ; "+names_all_func[real_s_ndx-1]+" - "+names_all_func[real_anchor_ndx-1]+" - "+names_all_func[real_H_ndx-1]+" m\n")
    angles_func.close()

def write_topology(top_file_func, bonds_file_func, angles_file_func):

    top_file = open(top_file_func, "r")
    cont_top_file = top_file.readlines()
    top_file.close()

    final_top=open(top_file_func, "w")
    final_top.close()
    final_top=open(top_file_func, "a")

    for i in range(len(cont_top_file)):
        final_top.writelines(cont_top_file[i])
        if ";   ai     aj funct   r             k" in cont_top_file[i]:
            bonds_file=open(bonds_file_func,"r")
            bonds_contents=bonds_file.readlines()
            bonds_file.close()
            for j in range(len(bonds_contents)):
                final_top.writelines(bonds_contents[j])
        if ";   ai     aj     ak    funct   theta         cth" in cont_top_file[i]:
            angles_file=open(angles_file_func,"r")
            angles_contents=angles_file.readlines()
            angles_file.close()
            for j in range(len(angles_contents)):
                final_top.writelines(angles_contents[j])

#Initializes names, types, coordinates, and number of atoms
names_all, types_all=init_topology(topfile_opt)
AU_xyz=init_element(grofile_opt, 'AU')
ST_xyz=init_element(grofile_opt, 'ST')
anchor_xyz=init_element(grofile_opt, anchorname_opt)
N_AU=len(AU_xyz)
N_ST=len(ST_xyz)
N_anchor=len(anchor_xyz)

#Calculates distance between all STs and all AUs, and STs and C1s
matrix_AU_ST=distance.cdist(ST_xyz, AU_xyz, 'euclidean')
matrix_anchor_ST=distance.cdist(ST_xyz, anchor_xyz, 'euclidean')

#Gets indexes of 2 closest AU atoms to each ST, also saving the unique indexes and number of times each unique index appears
all_mins=np.zeros((N_ST,2)).astype('int')
min_anchor=np.zeros(N_ST).astype('int')
for i in range(N_ST):
    all_mins[i,:]=matrix_AU_ST[i].argsort()[:2]
    min_anchor[i]=matrix_anchor_ST[i].argsort()[0]
unique_mins, counts_mins=np.unique(all_mins, return_counts=True)
unique_mins=unique_mins.astype('int')

#Gets coordinates of both close AUs, calculates distance, and angle
close_AU_xyz1=AU_xyz[all_mins[:,0],:]
close_AU_xyz2=AU_xyz[all_mins[:,1],:]
dist_AU_ST=calculate_distance(ST_xyz, close_AU_xyz1, close_AU_xyz2)
ang_AU_ST_AU=calculate_angles(ST_xyz, close_AU_xyz1, close_AU_xyz2)

#Makes list with as many Residue objects as STs
residues=[]
for i in range(N_ST):
    ST_res=np.array([i], dtype='int')
    AU_res=np.array(all_mins[i], dtype='int')
    for j in AU_res:
        common_ST=np.where(all_mins==j)[0]
        ST_res=np.append(ST_res, common_ST)
        AU_res=np.append(AU_res, all_mins[common_ST,:])
    n_au_bonds=np.array([], dtype='int')
    unique_AU_res=np.unique(AU_res)
    for j in unique_AU_res:
        n_au_bonds=np.append(n_au_bonds,counts_mins[unique_mins==j])
    residues.append(atom.Residue(number=i, type='WTF', gold_atoms=unique_AU_res, sulphur_atoms=np.unique(ST_res), gold_n_bonds=n_au_bonds))

#Residues with common atom lists are merged into one
new_residues=[]
res_number=1
AU_sets=[]
for i in range(len(residues)):
    res_act=residues[i]
    if not (np.any([(res_act.au_atoms == x).all() for x in AU_sets])):
        AU_sets.append(res_act.au_atoms)
        new_residues.append(atom.Residue(number=res_number, type='WTF', gold_atoms=res_act.au_atoms, sulphur_atoms=res_act.s_atoms, gold_n_bonds=res_act.au_n_bonds))
        res_number+=1

#Residues are assigned a type according to the AU-ST-AU angles present
for i in range(len(new_residues)):
    ns=len(new_residues[i].s_atoms)
    nau=len(new_residues[i].au_atoms)
    if ns==1 and nau==2:
        new_residues[i].change_type('STP')
    elif ns==2 and nau ==3:
        new_residues[i].change_type('STR')
    elif ns==3 and nau==4:
        angles_tmp=ang_AU_ST_AU[new_residues[i].s_atoms]
        if np.any(angles_tmp)<128.4 and np.any(angles_tmp) >110:
            new_residues[i].change_type('STV')
        elif np.any(angles_tmp)<=110 and np.any(angles_tmp) >90.8:
            new_residues[i].change_type('STC')

print_staple_pdb(new_residues, AU_xyz, ST_xyz, workdir_opt+'/staples.pdb')
print_new_bonds(new_residues, all_mins, min_anchor, workdir_opt+'/bonds.top', names_all, types_all, anchorname_opt)
print_new_angles(new_residues, all_mins, min_anchor, workdir_opt+'/angles.top', names_all, types_all, anchorname_opt, list_Hs_opt)
write_topology(topfile_opt, workdir_opt+'/bonds.top', workdir_opt+'/angles.top')
