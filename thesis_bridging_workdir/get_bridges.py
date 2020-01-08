"""
Written by: David S. Weber -

This script applies the Met-aromatic algorithm to a Protein Data Bank
entry of choosing, PDB_CODE. I wrote this script mainly for the completion
of my thesis to include in the code section. I feel that this is, to date,
my most elegant work. The script then finds bridging interactions.
"""

PDB_CODE = '1BPJ'  # input the PDB code of interest
D_C = 6.0  # this is the distance cutoff for the distance condition (2.3)
A_C = 109.5  # this is the angular cutoff for the angular condition (2.4)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# dependencies

import sys
import Bio.PDB as bp
import os
import pandas as pd
import numpy as np

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# constants and other hard coding 

STRUCT = 0
CHAIN = 'A'
ATOMS_MET = 'CE|SD|CG'  
ATOMS_TYR = 'CD1|CE1|CZ|CG|CD2|CE2'   
ATOMS_TRP = 'CD2|CE3|CZ2|CH2|CZ3|CE2' 
ATOMS_PHE = 'CD1|CE1|CZ|CG|CD2|CE2'   
LIST_ATOMS_TYR = ['CG', 'CD2', 'CE2', 'CZ', 'CE1', 'CD1']
LIST_ATOMS_TRP = ['CD2', 'CE3', 'CZ3', 'CH2', 'CZ2', 'CE2']
LIST_ATOMS_PHE = ['CG', 'CD2', 'CE2', 'CZ', 'CE1', 'CD1']
LIST_ATOMS_MET = ['CG', 'SD', 'CE']
TAGS = ['A', 'B', 'C', 'D', 'E', 'F']
SCAL1 = np.sin(np.pi / 2)
SCAL2 = 1 - np.cos(np.pi / 2)

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# function and class definitions

def vector_angle(u, v):
    # a routine for computing angle between two vectors
    num = np.dot(u, v)
    den = np.linalg.norm(v) * np.linalg.norm(u)
    return np.degrees(np.arccos(num / den))

def unit_vec(v):
    # get the unit vector of a vector
    return v / np.linalg.norm(v)
    
def angle_met_theta_rm(row):
    # df compatible routine for finding met-theta angle
    a = np.array([row['a_x_rm'], row['a_y_rm'], row['a_z_rm']])
    v = np.array([row['v_x'], row['v_y'], row['v_z']])
    return vector_angle(a, v)
    
def angle_met_phi_rm(row):
    # df compatible routine for finding met-phi angle
    g = np.array([row['g_x_rm'], row['g_y_rm'], row['g_z_rm']])
    v = np.array([row['v_x'], row['v_y'], row['v_z']])
    return vector_angle(g, v)
    
def angle_met_theta_cp(row):
    # df compatible routine for finding met-theta angle
    a = np.array([row['a_x_cp'], row['a_y_cp'], row['a_z_cp']])
    v = np.array([row['v_x'], row['v_y'], row['v_z']])
    return vector_angle(a, v)
    
def angle_met_phi_cp(row):
    # df compatible routine for finding met-phi angle
    g = np.array([row['g_x_cp'], row['g_y_cp'], row['g_z_cp']])
    v = np.array([row['v_x'], row['v_y'], row['v_z']])
    return vector_angle(g, v)

def res_segid(residue_input):
    # a workaround function for BioPython's get_segid() issue
    res_string = str(residue_input)
    resseq = res_string.find('resseq=') + 7
    return float(res_string[resseq:].split(' ')[0])
    
def get_hexagon_midpoints(x, y, z):
    """
    Function for computing midpoints between vertices in a hexagon
    Parameters:
        x, y, z -> list objects of x, y, and z hexagon coordinates
    Returns:
        x_mid, y_mid, z_mid -> a list of x, y, and z hexagon midpoint coordinates
    """
    
    # offset each list by 1
    x_f = x[1:] + [x[0]]
    y_f = y[1:] + [y[0]]
    z_f = z[1:] + [z[0]]

    # compute means between original and offset lists
    x_mid = [0.5 * (a + b) for a, b in zip(x, x_f)]
    y_mid = [0.5 * (a + b) for a, b in zip(y, y_f)]
    z_mid = [0.5 * (a + b) for a, b in zip(z, z_f)]
    
    return x_mid, y_mid, z_mid

class RodriguesMethod:
    """
    Here I use Rodrigues' rotation formula for completing the vertices
    of a regular tetrahedron. To start, we know two vertices of a tetrahedron,
    A and B, in addition to knowing the origin O. So we map A and B to the
    origin of the frame in which the tetrahedron resides by computing
    u = A - O and v = B - O. The directions of u, v are then flipped by scaling
    u, v by -1. We then take -u, -v and rotate the vectors by 90 degrees about
    k, where k is the line of intersection between the two orthogonal planes
    of the tetrahedron. These rotated vectors now describe the position of the
    remaining coordinates C, D. We have our tetrahedron with vertices A, B,
    C, D and the origin O.
    """
    def __init__(self, vertex_a, origin, vertex_b):
        # map to origin
        vec_u = vertex_a - origin
        vec_v = vertex_b - origin
        
        # flip direction
        self.vec_u = -1 * vec_u
        self.vec_v = -1 * vec_v
        
        # we then find the vector about which we rotate
        r = 0.5 * (self.vec_u + self.vec_v)
        
        # then find unit vector of r
        r_hat = r / np.linalg.norm(r)
        
        # get components of the unit vector r
        r_hat_x = r_hat[0]
        r_hat_y = r_hat[1]
        r_hat_z = r_hat[2]
        
        # get the W matrix
        W = np.matrix([[0, -r_hat_z, r_hat_y],
                       [r_hat_z, 0, -r_hat_x],
                       [-r_hat_y, r_hat_x, 0]])
                       
        # then construct Rodrigues rotation matrix
        self.R = np.matrix(np.eye(3)) + (SCAL1 * W) + (SCAL2 * np.matmul(W, W))
    
    # note that I flipped these methods to match previous algorithm
    def vector_g(self):
        # get vector g - first vertex
        return np.matmul(self.R, self.vec_u)

    def vector_a(self):
        # get vector v - second vertex
        return np.matmul(self.R, self.vec_v)

class LonePairs:
    # a class for computing the vectors parallel to MET SD lone pairs
    # methods associated with the class complete the vertices of a tetrahedron
    def __init__(self, terminal_a, midpoint, terminal_b):
        self.terminal_a = terminal_a
        self.midpoint = midpoint
        self.terminal_b = terminal_b
        self.u = self.terminal_a - self.midpoint # mapping to origin
        self.v = self.terminal_b - self.midpoint # mapping to origin
        self.NOT_vec = unit_vec(-0.5 * (unit_vec(self.v) + unit_vec(self.u)))
        
    def vector_a(self):
        cross_vec = np.cross(self.u, self.v)
        return self.NOT_vec + 2**0.5 * unit_vec(cross_vec)
        
    def vector_g(self):
        cross_vec = np.cross(self.v, self.u)
        return self.NOT_vec + 2**0.5 * unit_vec(cross_vec)
    

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# get data from Protein Data Bank

"""
We begin by importing the data of interest using BioPython. BioPython
has decent routines for connecting to the PDB and subsequently importing
the data of interest. I personally don't like BioPython for anything other
than importing and isolating PDB data so I use Pandas for all downstream work.
"""       

# create a biopython parser object
# QUIET=True suppresses 'WARNING: Chain [] is discontinuous at line...'
parser = bp.PDBParser(QUIET=True) 
pdbl = bp.PDBList()

# occasional issue with some PDB codes being decoded as numbers
try:
    path_to_file = pdbl.retrieve_pdb_file(pdb_code=PDB_CODE, file_format="pdb", pdir=os.getcwd())
except Exception as exception:
    sys.exit(exception)
else:
    structure = parser.get_structure('mystring', path_to_file)

# remove PDB file once we get structure
if os.path.exists(path_to_file):
    os.remove(path_to_file)     


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# apply the Met-aromatic algorithm in its entirety

"""
I also use BioPython to isolate PDB data down to the RA level (in SMCRA). 
The RA level is then further sorted using Pandas which is my 
preferred data analysis package.
"""                              

# isolate data of interest
model = structure[STRUCT]
chain = model[CHAIN]

# iteratively sort and arrange data for dataframe
to_df = []
for atoms in bp.Selection.unfold_entities(chain, 'A'):
    res_obj = atoms.get_parent()
    xyz = atoms.get_coord()
    to_df.append([chain.id, res_segid(res_obj), res_obj.get_resname(), atoms.name, xyz[0], xyz[1], xyz[2]])
                      
# cast to dataframe
df = pd.DataFrame(to_df)
df.columns = ['CHAIN', 'RES', 'RES NAME', 'ATOM', 'x', 'y', 'z']  

"""
This block is a vectorized equivalent of lines 1-5 in 
Algorithm 1 Met-aromatic. Here we isolate down to only those residues
of interest to us.
"""

# isolate down to residues we are interested in
df_MET = df[df['RES NAME'] == 'MET']
df_TYR = df[df['RES NAME'] == 'TYR']
df_TRP = df[df['RES NAME'] == 'TRP']  
df_PHE = df[df['RES NAME'] == 'PHE']     

"""
This block is a vectorized equivalent of lines 7-15 in
Algorithm 1 Met-aromatic. We strip down to only atoms of interest to us.
In this block we also order our data by tags in preparation for midpoint
calculations. This was not mentioned in Algorithm 1 Met-aromatic.
"""
           
# isolate down to atoms in the residues
df_MET = df_MET[df_MET['ATOM'].str.contains(ATOMS_MET)]
df_TYR = df_TYR[df_TYR['ATOM'].str.contains(ATOMS_TYR)]
df_TRP = df_TRP[df_TRP['ATOM'].str.contains(ATOMS_TRP)]
df_PHE = df_PHE[df_PHE['ATOM'].str.contains(ATOMS_PHE)]

# tag all atoms
df_MET['TAGS'] = df_MET['ATOM'].replace(LIST_ATOMS_MET, TAGS[0:3])
df_PHE['TAGS'] = df_PHE['ATOM'].replace(LIST_ATOMS_PHE, TAGS)
df_TYR['TAGS'] = df_TYR['ATOM'].replace(LIST_ATOMS_TYR, TAGS)
df_TRP['TAGS'] = df_TRP['ATOM'].replace(LIST_ATOMS_TRP, TAGS)

# sort atoms in preparation for midpoint calculation
df_MET = df_MET.sort_values(by='TAGS')
df_PHE = df_PHE.sort_values(by='TAGS')
df_TYR = df_TYR.sort_values(by='TAGS')
df_TRP = df_TRP.sort_values(by='TAGS')

# group by residue number
df_MET_GB = df_MET.groupby('RES')
df_PHE_GB = df_PHE.groupby('RES')
df_TYR_GB = df_TYR.groupby('RES')
df_TRP_GB = df_TRP.groupby('RES')

# first get all midpoints
midpoints_PHE = pd.DataFrame()
midpoints_TYR = pd.DataFrame()
midpoints_TRP = pd.DataFrame()

"""
Split / apply / combine procedure for computing hexagon midpoints. The df.append()
procedure, in general, is not recommended. However there are very few dfs to
append at this stage in the program and the df.append() function here has little
additional computational costs relative to other methods available.
"""

for residue_key, groupby_obj in df_PHE_GB:
    x_list = groupby_obj['x'].tolist()
    y_list = groupby_obj['y'].tolist()
    z_list = groupby_obj['z'].tolist()
    x_mid, y_mid, z_mid = get_hexagon_midpoints(x_list, y_list, z_list)
    dfs = pd.DataFrame([[residue_key] * 6, x_mid, y_mid, z_mid]).transpose()
    midpoints_PHE = midpoints_PHE.append(dfs)
    
for residue_key, groupby_obj in df_TYR_GB:
    x_list = groupby_obj['x'].tolist()
    y_list = groupby_obj['y'].tolist()
    z_list = groupby_obj['z'].tolist()
    x_mid, y_mid, z_mid = get_hexagon_midpoints(x_list, y_list, z_list)
    dfs = pd.DataFrame([[residue_key] * 6, x_mid, y_mid, z_mid]).transpose()
    midpoints_TYR = midpoints_TYR.append(dfs)

for residue_key, groupby_obj in df_TRP_GB:
    x_list = groupby_obj['x'].tolist()
    y_list = groupby_obj['y'].tolist()
    z_list = groupby_obj['z'].tolist()
    x_mid, y_mid, z_mid = get_hexagon_midpoints(x_list, y_list, z_list)
    dfs = pd.DataFrame([[residue_key] * 6, x_mid, y_mid, z_mid]).transpose()
    midpoints_TRP = midpoints_TRP.append(dfs)
    

# stack all midpoint dataframes
PHE_label = pd.DataFrame(midpoints_PHE.shape[0] * ['PHE'])
TYR_label = pd.DataFrame(midpoints_TYR.shape[0] * ['TYR'])
TRP_label = pd.DataFrame(midpoints_TRP.shape[0] * ['TRP'])

midpoints_PHE = midpoints_PHE.reset_index(drop=True)
midpoints_TYR = midpoints_TYR.reset_index(drop=True)
midpoints_TRP = midpoints_TRP.reset_index(drop=True)

midpoints_PHE = pd.concat([midpoints_PHE, PHE_label], axis=1)
midpoints_TYR = pd.concat([midpoints_TYR, TYR_label], axis=1)
midpoints_TRP = pd.concat([midpoints_TRP, TRP_label], axis=1)

midpoints_ALL = pd.concat([midpoints_PHE, midpoints_TYR, midpoints_TRP], axis=0)
midpoints_ALL = midpoints_ALL.reset_index(drop=True)
midpoints_ALL.columns = ['RES', 'x', 'y', 'z', 'ARO']


"""
Here we apply the lowest level of the Met-aromatic algorithm to our
refined feature space F. This corresponds to lines 17-26 in
Algorithm 1 Met-aromatic. Again, this is highly vectorized and does not
perfectly match the algorthm in the thesis.
"""

df_rm = pd.DataFrame()
df_cp = pd.DataFrame()

for res_key, met_res in df_MET.groupby('RES'):
    # first get CE, SD and CG coordinates
    CE, SD, CG = met_res.values[:, 4:7].astype(float)
    
    # then map all midpoints to origin of SD to get vectors v
    aro_res = pd.DataFrame()
    aro_res[['v_x', 'v_y', 'v_z']] = midpoints_ALL[['x', 'y', 'z']] - SD
    aro_res['RES'] = midpoints_ALL[['RES']].astype(int)
    aro_res['ARO'] = midpoints_ALL[['ARO']]
                 
    # then get norm of mapped vectors v
    aro_res['NORM'] = np.sqrt(np.square(aro_res[['v_x', 'v_y', 'v_z']]).sum(axis=1))
   
    # then apply distance condition & bank anything below cutoff
    aro_res = aro_res[aro_res['NORM'] <= D_C]
    
    # continue here to next iteration if aro_res is empty
    if aro_res.empty:
        continue
    else:
        pass
    
    # get lone pairs using rodrigues method approach
    obj_rotation = RodriguesMethod(CE, SD, CG)
    vec_a_rm = obj_rotation.vector_a().tolist()[0]
    vec_g_rm = obj_rotation.vector_g().tolist()[0]
    
    aro_res['a_x_rm'] = vec_a_rm[0]
    aro_res['a_y_rm'] = vec_a_rm[1]
    aro_res['a_z_rm'] = vec_a_rm[2]
    
    aro_res['g_x_rm'] = vec_g_rm[0]
    aro_res['g_y_rm'] = vec_g_rm[1]
    aro_res['g_z_rm'] = vec_g_rm[2]
    
    aro_res['MET-THETA-RM'] = aro_res.apply(angle_met_theta_rm, axis=1)
    aro_res['MET-PHI-RM'] = aro_res.apply(angle_met_phi_rm, axis=1)
    
    # get lone pairs using cross product method approach
    lp_object = LonePairs(CE, SD, CG)
    vec_a_cp = lp_object.vector_a()
    vec_g_cp = lp_object.vector_g()
    
    aro_res['a_x_cp'] = vec_a_cp[0]
    aro_res['a_y_cp'] = vec_a_cp[1]
    aro_res['a_z_cp'] = vec_a_cp[2]
    
    aro_res['g_x_cp'] = vec_g_cp[0]
    aro_res['g_y_cp'] = vec_g_cp[1]
    aro_res['g_z_cp'] = vec_g_cp[2]
    
    aro_res['MET-THETA-CP'] = aro_res.apply(angle_met_theta_cp, axis=1)
    aro_res['MET-PHI-CP'] = aro_res.apply(angle_met_phi_cp, axis=1)
    
    # now apply angular condition & bank anything meeting condition
    aro_res_rm = aro_res[(aro_res['MET-THETA-RM'] <= A_C) | (aro_res['MET-PHI-RM'] <= A_C)]
    aro_res_cp = aro_res[(aro_res['MET-THETA-CP'] <= A_C) | (aro_res['MET-PHI-CP'] <= A_C)]
    
    aro_res_rm = aro_res_rm[['ARO', 'RES', 'NORM', 'MET-THETA-RM', 'MET-PHI-RM']]
    aro_res_cp = aro_res_cp[['ARO', 'RES', 'NORM', 'MET-THETA-CP', 'MET-PHI-CP']]  
    
    aro_res_rm['MET'] = int(res_key)
    aro_res_cp['MET'] = int(res_key)
       
    # continue here to next iteration if aro_res_rm is empty
    if aro_res_rm.empty:
        continue
    else:
        df_rm = df_rm.append(aro_res_rm)
        
    # continue here to next iteration if aro_res_cp is empty
    if aro_res_cp.empty:
        continue
    else:
        df_cp = df_cp.append(aro_res_cp)
    
    
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# now find bridging interactions and output to console

df_rm = df_rm.reset_index(drop=True)
df_cp = df_cp.reset_index(drop=True)

print(df_rm)  # rodrigues method
print(df_cp)  # cross product method

print('-' * 50)
print('Bridges:')

rm_brdg = df_rm[['ARO', 'RES', 'MET']].drop_duplicates()
cp_brdg = df_cp[['ARO', 'RES', 'MET']].drop_duplicates()

rm_brdg['BRDG'] = rm_brdg.groupby('MET', as_index=False)['MET'].transform(lambda s: s.count())
cp_brdg['BRDG'] = cp_brdg.groupby('MET', as_index=False)['MET'].transform(lambda s: s.count())

print('Rodrigues method:')
print(rm_brdg)
print('Cross product method:')
print(cp_brdg)

