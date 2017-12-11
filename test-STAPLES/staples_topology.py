import numpy as np
from optparse import OptionParser
from scipy.spatial import distance
import atom
import matplotlib.pyplot as plt
parser=OptionParser()
parser.add_option("-x", "--gro", action="store", type='string', dest="GroFile", default='NP1.gro', help="Name of the .gro file with the staples")
parser.add_option("-p", "--top", action="store", type='string', dest="TopFile", default='NP1.top', help="Name of the .top file associated to the file in -x")
(options, args)= parser.parse_args()
grofile_opt=options.GroFile
topfile_opt=options.TopFile

def init_element(grofile_func, element):
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
    return np.array(elem_xyz_func), ini+1

def calculate_angles(pa, pb, pc):
    N_angles=len(pa[:,0])
    ab=np.subtract(pb, pa)
    ac=np.subtract(pc,pa)
    cos_ang=np.zeros(N_angles)
    for i in range(N_angles):
        cos_ang[i]=np.dot(ab[i,:], ac[i,:])/(np.linalg.norm(ab[i,:])*np.linalg.norm(ac[i,:]))
    return np.arccos(cos_ang)*180/np.pi

def calculate_distance(pa, pb, pc):
    N_dist=len(pa[:,0])
    dist_func=np.zeros((N_dist,2))
    pab=np.subtract(pb,pa)
    pac=np.subtract(pc,pa)
    for i in range(N_dist):
        dist_func[i,0]=np.linalg.norm(pab[i,:])
        dist_func[i,1]=np.linalg.norm(pac[i,:])
    return dist_func

def print_xyz(coordenadas, nombres, fnombre):
    output=open(fnombre, "w")
    output.write(str(len(nombres)) + "\n \n")
    for i in range(len(nombres)):
        output.write(str(nombres[i]) + " " + str(coordenadas[i,0]) + " " + str(coordenadas[i,1]) + " " + str(coordenadas[i,2]) + "\n")
    output.close()

AU_xyz, AU_ini_at=init_element(grofile_opt, 'AU')
ST_xyz, ST_ini_at=init_element(grofile_opt, 'ST')
N_AU=len(AU_xyz)
N_ST=len(ST_xyz)

matrix_AU_ST=distance.cdist(ST_xyz, AU_xyz, 'euclidean')

all_mins=np.zeros((N_ST,2))
for i in range(N_ST):
    mins_ndx=matrix_AU_ST[i].argsort()[:2]
    all_mins[i,:]=mins_ndx
all_mins=all_mins.astype('int')
unique_mins, counts_mins=np.unique(all_mins, return_counts=True)
unique_mins=unique_mins.astype('int')

close_AU_xyz1=AU_xyz[all_mins[:,0],:]
close_AU_xyz2=AU_xyz[all_mins[:,1],:]
dist_AU_ST=calculate_distance(ST_xyz, close_AU_xyz1, close_AU_xyz2)
ang_AU_ST_AU=calculate_angles(ST_xyz, close_AU_xyz1, close_AU_xyz2)

ST_check=np.zeros(N_ST)
for i in range(N_ST):
    if ST_check[i]==0:
        ST_check[i]=1
        ST_res=[]
        AU_res=[]

        ST_res.append(i)
        AU_res.append(all_mins[i,:])
