# written by dsw7@sfu.ca
# automates entire DFT procedure for H2S / C6H6 angular / distance translation 

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib           import cm
from xyz_file_generator   import get_XYZ_file
from numpy                import array, matrix, linspace, radians, cos, sin
from os                   import path, listdir, mkdir, chdir, stat, chmod
from re                   import search
from subprocess           import call
from stat                 import S_IEXEC


# make a folder to dump all .xyz data
if not path.exists('orca_workdir'):
    mkdir('orca_workdir')

def rotY(a, x, y, z):
    # y axis rotation + translation homogeneous transform
    return matrix([
    [ cos(a), 0,     sin(a), x],
    [         0, 1,       0, y],
    [-sin(a), 0,     cos(a), z],
    [         0, 0,       0,  1]])


# *******************************************************************
# prepare benzene coordinates

coords_C6H6 = [
  ['C',        0.00000,        1.40272,        0.00000],
  ['H',        0.00000,        2.49029,        0.00000],
  ['C',       -1.21479,        0.70136,        0.00000],
  ['H',       -2.15666,        1.24515,        0.00000],
  ['C',       -1.21479,       -0.70136,        0.00000],
  ['H',       -2.15666,       -1.24515,        0.00000],
  ['C',        0.00000,       -1.40272,        0.00000],
  ['H',        0.00000,       -2.49029,        0.00000],
  ['C',        1.21479,       -0.70136,        0.00000],
  ['H',        2.15666,       -1.24515,        0.00000],
  ['C',        1.21479,        0.70136,        0.00000],
  ['H',        2.15666,        1.24515,        0.00000]
]


# *******************************************************************
# prepare hydrogen sulfide coordinates

ANGLE_TETRAHEDRON = 109.5
CORRECTION = ANGLE_TETRAHEDRON - 92.1 
LP_ANGLE = ANGLE_TETRAHEDRON + CORRECTION
HALF_LP_ANGLE = 0.5 * LP_ANGLE

# raw H2S coordinates
p1 = array([[0.927], [ 0.962], [0.000], [1.000]])
p2 = array([[0.000], [ 0.000], [0.000], [1.000]])
p3 = array([[0.927], [-0.962], [0.000], [1.000]])


# *******************************************************************
# rotate + translate H2S coordinates

tr_start = 3.50   # where to start translation (Angstroms)
tr_end   = 4.50  # where to end translation (Angstroms)
tr_incre = 0.50  # increment between translations

theta_start = 0.00   # where to start angle (degrees)
theta_end   = 45.0   # where to end angle (degrees)
theta_incre = 5.00   # increment

for t in linspace(tr_start, tr_end, int(1 + (tr_end - tr_start) / tr_incre)):
    if t - 1.21479 > 6.00:
        break  # if ||v|| exceeds 6.0 angstroms
    else:
        for theta in linspace(theta_start, 
                              theta_end, 
                              int(1 + (theta_end - theta_start) / theta_incre)):
            # perform the homogeneous transforms
            rot_p1 = rotY(radians(HALF_LP_ANGLE + theta), t, 0.00, 0.00) * p1
            rot_p2 = rotY(radians(HALF_LP_ANGLE + theta), t, 0.00, 0.00) * p2
            rot_p3 = rotY(radians(HALF_LP_ANGLE + theta), t, 0.00, 0.00) * p3
            
            # prepare for export
            rot_p1 = ['H'] + rot_p1[0:3].transpose().tolist()[0]
            rot_p2 = ['S'] + rot_p2[0:3].transpose().tolist()[0]
            rot_p3 = ['H'] + rot_p3[0:3].transpose().tolist()[0]
            
            coords_all = coords_C6H6 + [rot_p1, rot_p2, rot_p3]
            get_XYZ_file(coords_all, 'orca_workdir/all_{}_{}_.xyz'.format(t - 1.21479, theta))
        

# *******************************************************************
# make .inp decks        
        
# scope out pwd to find .xyz files
local_dirs = listdir('orca_workdir')
local_dirs = [i for i in local_dirs if path.splitext(i)[1] == '.xyz']

# using a docstring fails
HEADER = [
'## avogadro generated ORCA input file\n',
'# Advanced Mode\n',
'#\n'
'! BP RI SP def2-SVP NormalPrint Grid4  NormalSCF\n',
'%scf\n'
'	MaxIter 125\n',
'	CNVDIIS 1\n',
'	CNVSOSCF 1\n',
'end\n',
'%output\n',
'	print[p_mos] true\n',
'	print[p_basis] 5\n',
'end\n',

'* xyz 0 1\n']

FOOTER = ['*\n']

for file in local_dirs:
    # read in all .xyz files
    with open(path.join('orca_workdir', file), 'r') as f_in:
        lines = f_in.readlines()

    lines = [item for item in lines if search(r'C|S|H', item)]  # strip .xyz header

    input_filename = file.split('.')
    output_filename = '{}.{}.inp'.format(*input_filename)
    
    # export all .inp files
    with open(path.join('orca_workdir', output_filename), 'w') as f_out:
        for line in HEADER + lines + FOOTER:
            f_out.write(line)

        
# *******************************************************************
# make shell file
    
# scope out pwd to find .xyz files
local_dirs = listdir('orca_workdir')

# generate input and output deck names
files_input = [i for i in local_dirs if path.splitext(i)[1] == '.inp']
files_output = [i.replace('.inp', '.out') for i in files_input]

in_str = ''
for f in files_input:
    in_str += f + '\n'

array_input = """inp_decks=({}) \n""".format(in_str)

out_str = ''
for f in files_output:
    out_str += f + '\n'

array_output = """out_decks=({}) \n""".format(out_str)

loop = """
for((i=0; i<={}; i++));
do
    orca ${{inp_decks[i]}} > ${{out_decks[i]}}
done
""".format(len(files_input))

# merge all code
output_sh_code = """#!/bin/sh\n""" + array_input + array_output + loop

# write to a shell file
with open(r'orca_workdir/inp_to_out.sh', 'w') as shellfile:
    shellfile.write(output_sh_code)


# *******************************************************************
# execute shell file  

chdir('orca_workdir')
st = stat('inp_to_out.sh')
chmod('inp_to_out.sh', st.st_mode | S_IEXEC)
rc = call([r'./inp_to_out.sh'])

# appears to be an issue with termination here
# orca continues processing despite closing python script
# query the status of the shell command?



# *******************************************************************
# read in all energy levels

outdirs = listdir()
outdirs = [i for i in outdirs if path.splitext(i)[1] == '.out']

HOMO_energies = []
LUMO_energies = []

# iterate over all .out files
for file in outdirs:
    with open(file) as f:
        lines = f.read().split('\n')
    
    # get index of where orbital energies start
    idx =  lines.index('ORBITAL ENERGIES')
    
    # get the first 100 lines
    energies = lines[idx + 4: idx + 104]
    
    # split by whitespace
    energies = [e.split() for e in energies]
    
    # we want indices 1 and 3 -> 0CC & E (eV)
    CC = [i[1] for i in energies]
    eV = [i[3] for i in energies]
    
    # we want to get the last index of 2.00 in CC - this is where HOMO -> LUMO
    idx_switch = 100 - CC[::-1].index('2.0000')
    
    """
    CC[idx_switch - 1])  # the last 2.0000
    CC[idx_switch])      # the first 0.0000
    CC[idx_switch + 1])  # the second 0.0000
    """
    
    # get HOMO + LUMO
    HOMO = float(eV[idx_switch - 1])
    LUMO = float(eV[idx_switch])
    
    fileparts = file.split('_')
    
    # get ||v|| from filename
    norm_v = float(fileparts[1])
    
    # get angle
    angle = float(fileparts[2].split('.')[0])
    
    HOMO_energies.append((norm_v, angle, HOMO))
    LUMO_energies.append((norm_v, angle, LUMO))


# sort first by ||v|| then by incremented angle
HOMO_energies = sorted(HOMO_energies, key=lambda a: a[0:2])
LUMO_energies = sorted(LUMO_energies, key=lambda a: a[0:2])  

# transpose to isolate x, y, z data
HOMO_transposed = array(HOMO_energies).transpose()
LUMO_transposed = array(LUMO_energies).transpose()

# *******************************************************************
# the final plot

fig = plt.figure()
ax = fig.gca(projection='3d')
ax.plot_trisurf(*HOMO_transposed, cmap=cm.jet, linewidth=0.2)
