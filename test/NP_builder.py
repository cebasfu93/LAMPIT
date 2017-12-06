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
outname_opt=options.OutputFile
corename_opt=options.CoreFile
ligname_opt=options.LigandFile
res_name_opt=options.ResidueName

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
d_CS=1.83 #Angstroms (wikipedia)

ar_const=np.zeros(N_S)
ar_mag=np.zeros(N_S)
ar_s=np.zeros((N_S,3))
ar_new=np.zeros((N_S,3))
def center_COM(xyz_obj_func, names_core_func, objeto_COM):
    COM=np.average(xyz_obj_func[names_core_func==objeto_COM,:], axis=0)
    for i in range(len(xyz_obj_func[:,0])):
        xyz_obj_func[i,:]=xyz_obj_func[i,:]-COM
    return xyz_obj_func

def get_new_linker(xyz_S_func, dist):
    mag_S=np.linalg.norm(xyz_S_func)
    const=(mag_S+dist)/mag_S
    new_linker_tmp=const*xyz_S_func
    return const, mag_S, xyz_S_func, new_linker_tmp

def place_linker(xyz_core_func, names_core_func, xyz_lig_func, terminal_ndx, objeto_link):
    xyz_S_func=xyz_core_func[np.where(names_core_func==objeto_link)[0],:]
    new_linkers_func=np.zeros((N_S, 3))
    old_linker_func=xyz_lig_func[terminal_ndx,:]
    for i in range(N_S):
        ar_const[i], ar_mag[i], ar_s[i], ar_new[i]=get_new_linker(xyz_S_func[i,:], d_CS)
        new_linker_func=ar_new[i]
        new_linkers_func[i]=new_linker_func

    """fig=plt.figure()
    plt.subplot(2,2,1)
    plt.title('constants')
    plt.plot(ar_const)
    plt.subplot(2,2,2)
    plt.title('mag_S')
    plt.plot(ar_mag)
    plt.subplot(2,2,3)
    plt.title('s position')
    plt.plot(ar_s)
    plt.subplot(2,2,4)
    plt.title('new_C')
    plt.plot(ar_new)
    plt.show()"""
    return new_linkers_func

def place_coating(linked_NP_func, xyz_lig_func):
    N_lig_tot=N_at_lig*N_S
    coating_func=np.zeros((N_lig_tot,3))
    for i in range(N_S):
        coating_func[i*N_at_lig]=linked_NP_func[i]
        for j in range(N_at_lig):
            coating_func[i*N_at_lig+j]=xyz_lig_func[j]
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
linkers=place_linker(centered_core, names_core, xyz_lig, 0, 'S')
ligands=place_coating(linkers, xyz_lig)
print_NP_pdb(names_lig, names_core, ligands, centered_core, types_lig, outname_opt)
