#imports necessary packages and classes
import numpy as np
from optparse import OptionParser
from scipy.spatial import distance
import atom
import matplotlib.pyplot as plt

#Declaration of the flags to run the script
parser=OptionParser()
parser.add_option("-x", "--gro", action="store", type='string', dest="GroFile", default='NP1.gro', help="Name of the .gro file with the staples")
parser.add_option("-p", "--top", action="store", type='string', dest="TopFile", default='NP1.top', help="Name of the .top file associated to the file in -x")
(options, args)= parser.parse_args()
grofile_opt=options.GroFile
topfile_opt=options.TopFile

def init_element(grofile_func, element):
    #Returns the xyz coordinates of the atoms with name 'element'
    elem_file_func=np.genfromtxt(grofile_func, delimiter='\n', dtype='str', skip_header=2)
    first=True
    ini=0
    fin=0
    N_elem_func=0
    for i in range(len(elem_file_func)):
        if element in elem_file_func[i] and first:
            ini=i
            first=False
        if element not in elem_file_func[i] and not first:
            fin=i
            break
    N_elem_func=fin-ini
    elem_xyz_func=np.zeros((N_elem_func,3))
    elem_names_func=[]
    for i in range(N_elem_func):
        elem_actual=elem_file_func[ini+i].split()
        elem_xyz_func[i,:]=elem_actual[4:7]
        elem_names_func.append(elem_actual[2])
    return np.array(elem_xyz_func)*10, ini+1

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
        types_all_func.append(types_file_func[i].split()[1])
    return names_all_func, types_all_func

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
        for j in range(len(res_act.s_atoms)):
            at+=1
            write_pdb_block('ST', res_act.restype, ST_xyz_func[res_act.s_atoms[j],:], res_act.resnum, at, out_fname)
        for j in range(len(res_act.au_atoms)):
            at+=1
            write_pdb_block('AU', res_act.restype, AU_xyz_func[res_act.au_atoms[j],:], res_act.resnum, at, out_fname)
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

def print_new_bonds(res_list_tmp, AU_ini_at_tmp, ST_ini_at_tmp):

    return 0
    #Prints the staple-staple and staple-ligand bonds
    return 0
def print_new_angles():
    #Prints the angles with staple atoms present
    return 0

def print_new_dihedrals():
    return 0
    #Prints the proper dihedrals with staple atoms present
    return 0

def write_topology(topfile_func, res_list, AU_ini_at_func, ST_ini_at_func, names_all_func):
    #Prints the complete topology of the system with the new parameters
    top_file_func=np.genfromtxt(topfile_func, delimiter='\n', dtype='str')
    for i in range(len(top_file_func)):
        #print(top_file_func[i])
        if ";   ai     aj funct   r             k" in top_file_func[i]:
            print_new_bonds(res_list, AU_ini_at_func, ST_ini_at_func)
        if ";   ai     aj     ak    funct   theta         cth" in top_file_func[i]:
            print_new_angles()
            print(i)
        if ";    i      j      k      l   func   phase     kd      pn" in top_file_func[i]:
            print_new_dihedrals()
            print(i)

#Initializes names, types, coordinates, and number of atoms
names_all, types_all=init_topology(topfile_opt)
AU_xyz, AU_ini_at=init_element(grofile_opt, 'AU')
ST_xyz, ST_ini_at=init_element(grofile_opt, 'ST')
N_AU=len(AU_xyz)
N_ST=len(ST_xyz)

#Calculates distance between all STs and all AUs
matrix_AU_ST=distance.cdist(ST_xyz, AU_xyz, 'euclidean')

#Gets indexes of 2 closest AU atoms to each ST, also saving the unique indexes and number of times each unique index appears
all_mins=np.zeros((N_ST,2))
for i in range(N_ST):
    mins_ndx=matrix_AU_ST[i].argsort()[:2]
    all_mins[i,:]=mins_ndx
all_mins=all_mins.astype('int')
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
    residues.append(atom.Residue(number=i, type='WTF', gold_atoms=np.unique(AU_res), sulphur_atoms=np.unique(ST_res)))
new_residues=[]
res_number=1
AU_sets=[]
#Residues with common atom lists are merged into one
for i in range(len(residues)):
    if not (np.any([(residues[i].au_atoms == x).all() for x in AU_sets])):
        AU_sets.append(residues[i].au_atoms)
        new_residues.append(atom.Residue(number=res_number, type='WTF', gold_atoms=residues[i].au_atoms, sulphur_atoms=residues[i].s_atoms))
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
            new_residues[i].change_type('STC')
        elif np.any(angles_tmp)<=110 and np.any(angles_tmp) >90.8:
            new_residues[i].change_type('STV')

#print_staple_pdb(new_residues, AU_xyz, ST_xyz, 'test.pdb')
write_topology(topfile_opt, new_residues, AU_xyz, ST_xyz, names_all)
