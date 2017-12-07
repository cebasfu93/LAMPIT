import numpy as np
import matplotlib.pyplot as plt
from optparse import OptionParser
from transformations import affine_matrix_from_points, superimposition_matrix


parser=OptionParser()
parser.add_option("-o", "--output", action="store", type='string', dest="OutputFile", default='Au144L1-60.xyz', help="Name of the output file")
parser.add_option("-c", "--core", action="store", type='string', dest="CoreFile", default='Au144_MD6341surface1.xyz', help="Name of the xyz core file (includes S)")
parser.add_option("-l", "--ligand", action="store", type='string', dest="LigandFile", default='ligand.mol2', help="Name of the mol2 uncapped ligand file")
parser.add_option("-r", "--resname", action="store", type='string', dest="ResidueName", default='UNK', help="3-letter name of the coating ligand")
parser.add_option("-a", "--anchor", action="append", dest="AnchorIndexes", default=[], help="Indexes (starting from 0) of the atoms that will be aligned to the COM-S vector")
(options, args)= parser.parse_args()
outname_opt=options.OutputFile
corename_opt=options.CoreFile
ligname_opt=options.LigandFile
res_name_opt=options.ResidueName
anchor_ndx=np.array(list(map(int, options.AnchorIndexes[0].split(','))))

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
    names_lig_func=[]
    types_lig_func=[]
    xyz_lig_func=np.zeros((n_at,3))
    for i in range(n_at):
        at_file=lig_file[i].split()
        names_lig_func.append(at_file[1])
        types_lig_func.append(at_file[5])
        xyz_lig_func[i,:]=at_file[2:5]
    return np.array(xyz_lig_func), np.array(names_lig_func), np.array(types_lig_func)

def init_core_pdb(core_fname):
    core_file=np.genfromtxt(core_fname, delimiter='\n', dtype=str)
    first=1
    N_core_func=0
    for i in range(len(core_file)):
        if "ATOM" in core_file[i]:
            N_core_func+=1
            if first==1:
                ini=i
                first=2
    core_file=core_file[ini:ini+N_core_func]
    names_core_func=[]
    xyz_core_func=np.zeros((N_core_func,3))
    for i in range(N_core_func):
        at_file=core_file[i].split()
        names_core_func.append(at_file[2])
        xyz_core_func[i,:]=at_file[5:8]
    return np.array(xyz_core_func), np.array(names_core_func)

xyz_lig, names_lig, types_lig=init_lig_mol2(ligname_opt)
xyz_core, names_core=init_core_pdb(corename_opt)

N_at_lig=len(xyz_lig[:,0])
N_S=len(names_core[names_core=='S'])
N_anchors=len(anchor_ndx)
d_CS=1.83 #Angstroms (wikipedia)

def center_COM(xyz_obj_func, names_core_func, objeto_COM):
    COM=np.average(xyz_obj_func[names_core_func==objeto_COM,:], axis=0)
    for i in range(len(xyz_obj_func[:,0])):
        xyz_obj_func[i,:]=xyz_obj_func[i,:]-COM
    return xyz_obj_func

def get_linkers_C(xyz_core_func, names_core_func, name_linker_tmp):
    return xyz_core_func[names_core_func==name_linker_tmp,:]

def get_stone(xyz_stones_tmp, xyz_lig_tmp, stone_ndx):
    mag_C=np.linalg.norm(xyz_stones_tmp)
    scaling=(mag_C+np.linalg.norm(xyz_lig_tmp[anchor_ndx[stone_ndx],:]-xyz_lig_tmp[anchor_ndx[0],:]))/mag_C
    return scaling*xyz_stones_tmp

def place_stones(xyz_core_func, names_core_func, xyz_lig_func, names_lig_func, name_linker_func):
    xyz_link_C=get_linkers_C(xyz_core_func, names_core_func, name_linker_func)
    N_stones_func=N_anchors*N_S
    xyz_stones_func=np.zeros((N_stones_func,3))
    names_stones_func=np.zeros(N_stones_func)
    names_stones_func=names_stones_func.astype('str')
    for i in range(N_S):
        xyz_stones_func[i*N_anchors,:]=xyz_link_C[i,:]
        names_stones_func[i*N_anchors]=names_lig_func[anchor_ndx[0]]
        for j in range(1,N_anchors):
            xyz_stones_func[i*N_anchors+j,:]=get_stone(xyz_stones_func[i*N_anchors,:], xyz_lig_func, j)
            names_stones_func[i*N_anchors+j]=names_lig_func[anchor_ndx[j]]
    return xyz_stones_func, names_stones_func

def coat_NP(xyz_core_func, names_core_func, xyz_stones_func, names_stones_func):
    keep_rows=[]
    for i in range(len(names_core_func)):
        if names_core_func[i]=='Au' or names_core_func[i]=='S':
            keep_rows.append(i)
    xyz_naked_func=xyz_core_func[keep_rows,:]
    names_naked_func=names_core_func[keep_rows]

    xyz_coated_func=np.append(xyz_naked_func, xyz_stones_func, axis=0)
    names_coated_func=np.append(names_naked_func, names_stones_func, axis=0)
    return xyz_coated_func, names_coated_func

def print_NP_pdb(xyz_coated_func, names_coated_func, names_lig_func, out_fname):
    output=open(out_fname, "w")
    res=1
    at=0
    for i in range(len(names_coated_func)):
        at+=1
        if names_coated_func[i]=='Au' or names_coated_func[i]=='S':
            write_pdb_block(names_coated_func[i], 'COR', xyz_coated_func[i,:], res, at, out_fname)
        else:
            if(names_coated_func[i]==names_lig_func[anchor_ndx[0]]):
                res+=1
            write_pdb_block(names_coated_func[i], res_name_opt, xyz_coated_func[i,:], res, at, out_fname)
    output.close()

def write_pdb_block(atname_func, res_name_func, xyz_func, resnum, atnum, out_filename):
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
    output=open(fnombre, "w")
    output.write(str(len(nombres)) + "\n \n")
    for i in range(len(nombres)):
        output.write(str(nombres[i]) + " " + str(coordenadas[i,0]) + " " + str(coordenadas[i,1]) + " " + str(coordenadas[i,2]) + "\n")
    output.close()

centered_core=center_COM(xyz_core, names_core, 'Au')
xyz_stones, names_stones=place_stones(centered_core, names_core, xyz_lig, names_lig, 'C')
xyz_coated_NP, names_coated_NP=coat_NP(centered_core, names_core, xyz_stones, names_stones)
#print_xyz(xyz_coated_NP, names_coated_NP, outname_opt)
print_NP_pdb(xyz_coated_NP, names_coated_NP, names_lig, outname_opt)
