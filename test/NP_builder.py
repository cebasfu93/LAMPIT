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
    return xyz_lig_func, names_lig_func, types_lig_func

def init_core_xyz(core_fname):
    core_file=np.genfromtxt(core_fname, skip_header=2, dtype='str')
    names_core_func=core_file[:,0]
    xyz_core_func=core_file[:,1:].astype('float')
    return xyz_core_func, names_core_func

xyz_lig, names_lig, types_lig=init_lig_mol2(ligname_opt)
xyz_core, names_core=init_core_xyz(corename_opt)

N_at_lig=len(xyz_lig[:,0])
N_S=len(names_core[names_core=='S'])
N_anchors=len(anchor_ndx)
d_CS=1.83 #Angstroms (wikipedia)

def center_COM(xyz_obj_func, names_core_func, objeto_COM):
    COM=np.average(xyz_obj_func[names_core_func==objeto_COM,:], axis=0)
    for i in range(len(xyz_obj_func[:,0])):
        xyz_obj_func[i,:]=xyz_obj_func[i,:]-COM
    return xyz_obj_func

def get_new_anchors(xyz_S_tmp, xyz_lig_tmp, dist, anchor_ndx_tmp):
    new_anchors_tmp=np.zeros((N_anchors,3))
    for i in range (N_anchors):
        mag_S=np.linalg.norm(xyz_S_tmp)
        if i==0:
            const=(mag_S+dist)/mag_S
            new_anchors_tmp[i,:]=const*xyz_S_tmp
        else:
            rel_dist=np.linalg.norm(xyz_lig_tmp[anchor_ndx_tmp[i],:]-xyz_lig_tmp[anchor_ndx_tmp[i-1],:])
            const=(mag_S+dist+rel_dist)/mag_S
            new_anchors_tmp[i,:]=const*xyz_S_tmp
    return new_anchors_tmp

def place_anchors(xyz_core_func, names_core_func, xyz_lig_func, objeto_link, anchor_ndx_func):
    xyz_S_func=xyz_core_func[np.where(names_core_func==objeto_link)[0],:]
    N_tot_anchors=N_anchors*N_S
    new_anchors_func=np.zeros((N_tot_anchors, 3))
    old_anchors_func=xyz_lig_func[anchor_ndx_func,:]
    for i in range(N_S):
        new_anchors_func[i*N_anchors:(i+1)*N_anchors,:]=get_new_anchors(xyz_S_func[i,:], xyz_lig_func, d_CS, anchor_ndx_func)
    return old_anchors_func, new_anchors_func

def place_coating(new_anchors_func, xyz_lig_func, anchor_ndx_func):
    old_anchors_func=xyz_lig_func[anchor_ndx_func,:]
    N_lig_tot=N_at_lig*N_S
    coating_func=np.zeros((N_lig_tot,3))
    for i in range(N_S):
        trans_matrix=affine_matrix_from_points(old_anchors_func, new_anchors_func[i*N_anchors:(i+1)*N_anchors,:], shear=False, scale=False)
        for j in range(N_at_lig):
            if j in anchor_ndx_func:
                coating_func[i*N_at_lig+j]=new_anchors_func[i*N_anchors+np.where(anchor_ndx_func==j)[0][0],:]
            else:
                victime_vec=np.append(xyz_lig_func[j,:],1)
                coating_func[i*N_at_lig+j]=np.dot(trans_matrix, victime_vec)[0:3]
    return coating_func


def print_NP_pdb(names_lig_func, names_core_func, xyz_ligs_func, xyz_core_func, types_lig_func, out_fname):
    output=open(out_fname, "w")
    res=1
    at=0
    for i in range(len(xyz_core_func[:,0])):
        at+=1
        write_pdb_block(names_core_func[i], 'COR', xyz_core_func[i,:], res, at, out_fname)
    for i in range(N_S):
        res+=1
        for j in range(N_at_lig):
            at+=1
            write_pdb_block(names_lig_func[j], res_name_opt, xyz_ligs_func[i*N_at_lig+j,:], res, at, out_fname)
    output.close()

def write_pdb_block(atname_func, res_name_func, xyz_func, resnum, atnum, out_filename):
    xyz_func=np.round(xyz_func, decimals=4)
    coords=open(out_filename, 'a')
    coords.write('ATOM'.ljust(6))
    coords.write(str(atnum).rjust(5))
    coords.write('  '+str(atname_func).ljust(3))
    coords.write(' '+str(res_name_func).ljust(3))
    coords.write('  '+str(resnum).rjust(4))
    coords.write('   '+ str(xyz_func[0]).rjust(8))
    coords.write(str(xyz_func[1]).rjust(8))
    coords.write(str(xyz_func[2]).rjust(8)+"\n")
    coords.close()

centered_core=center_COM(xyz_core, names_core, 'Au')
old_anchors, new_anchors=place_anchors(centered_core, names_core, xyz_lig, 'S', anchor_ndx)
ligands=place_coating(new_anchors, xyz_lig, anchor_ndx)
print_NP_pdb(names_lig, names_core, ligands, centered_core, types_lig, outname_opt)
