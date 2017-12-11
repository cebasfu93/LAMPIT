import numpy as np
from optparse import OptionParser
from scipy.spatial import distance

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
    return elem_xyz_func, ini+1

AU_xyz, AU_ini_at=init_element(grofile_opt, 'AU')
ST_xyz, ST_ini_at=init_element(grofile_opt, 'ST')

dist_AU_ST=distance.cdist(ST_xyz, AU_xyz, 'euclidean')
