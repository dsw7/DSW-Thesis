"""
This script fetches a Met-aromatic interaction of choosing from the PDB
and then generates a set of translations for the aromatic interaction.
"""

PDB_CODE = '1cfc'  # input the PDB code of interest
AROMATIC_RESIDUE = 12
METHIONINE_RESIDUE = 72
TRANS_START = -1.0
TRANS_END = 0.75
INCREMENT = 9


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

def res_segid(residue_input):
    # a workaround function for BioPython's get_segid() issue
    res_string = str(residue_input)
    resseq = res_string.find('resseq=') + 7
    return int(res_string[resseq:].split(' ')[0])
    
def get_XYZ_file(input_list, file_name):
    # I inherited this function from the now defunct PDBFileDriver library which
    # I wrote in early 2017 - still has many useful functions
    """
    PDBFileDriver.get_XYZ_file(input_list)
    
    A base script for generating .xyz Avogadro input files.
    
    Parameters : input_list, a list of elements and 3d coordinates formatted as follows:
                 input_list = [
                               ['N',1.00,2.00,1.00],
                               ['H',1.00,4.00,1.00], ...
                              ]
                 file_name, name of export file.
    
    The script generates:
        xyz_data.xyz
    """
    file = open(file_name + '.xyz', 'w')   
    file.write(str(len(input_list)) + '\n')
    file.write('PDBFileDriver generated XYZ file\n')    
    for i in range(0, len(input_list)):
        string = ' {} {} {} {} \n'.format(*i)
        file.write(string)
    file.close()

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

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
                         

# isolate data of interest
model = structure[STRUCT]
chain = model[CHAIN]

# iteratively sort and arrange data for dataframe
to_df = []
for atoms in bp.Selection.unfold_entities(chain, CHAIN):
    res_obj = atoms.get_parent()
    xyz = atoms.get_coord()
    to_df.append([chain.id, res_segid(res_obj), res_obj.get_resname(), atoms.name, xyz[0], xyz[1], xyz[2]])
          
                   
# cast to dataframe
df = pd.DataFrame(to_df)
df.columns = ['CHAIN', 'RES', 'RES NAME', 'ATOM', 'x', 'y', 'z']  


# isolate down to residue
df_MET = df[df['RES'] == METHIONINE_RESIDUE]
df_ARO = df[df['RES'] == AROMATIC_RESIDUE]


# isolate down to atoms  
met_ID = df_MET['RES NAME'].drop_duplicates().iloc[0]
aro_ID = df_ARO['RES NAME'].drop_duplicates().iloc[0]


# methionine atoms
if met_ID == 'MET':
    df_MET = df_MET[df_MET['ATOM'].str.contains(ATOMS_MET)]
else:
    sys.exit('Residue {} is not MET.'.format(METHIONINE_RESIDUE))
    

# aromatic atoms
if aro_ID == 'PHE':
    df_ARO = df_ARO[df_ARO['ATOM'].str.contains(ATOMS_PHE)]
elif aro_ID == 'TYR':
    df_ARO = df_ARO[df_ARO['ATOM'].str.contains(ATOMS_TYR)]
elif aro_ID == 'TRP':
    df_ARO = df_ARO[df_ARO['ATOM'].str.contains(ATOMS_TRP)]
else:
    sys.exit('Residue {} is not any of PHE, TYR or TRP.'.format(AROMATIC_RESIDUE))
    
    
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# we have our raw data - now perform translation
    
vectors = df_ARO[['x', 'y', 'z']].values  # vectors parallel to aromatic plane

u = vectors[2] - vectors[0]
v = vectors[1] - vectors[0]

n = np.cross(u, v)  # the vector normal to the aromatic plane
n_hat = n / np.linalg.norm(n)  # get unit vector of normal vector

df_MET['w'] = 1.0  # will need this downstream

for INCR_TRANS in np.linspace(TRANS_START, TRANS_END, INCREMENT):
    n_hat_scaled = INCR_TRANS * n_hat  # scale unit vector to translation increment
    
    # place n_hat_scaled vector coordinates into homogeneous transformation matrix
    x_t = n_hat_scaled[0]
    y_t = n_hat_scaled[1]
    z_t = n_hat_scaled[2]
    
    # our transformation matrix...
    trans_matrix = np.matrix([[1, 0, 0, x_t],
                              [0, 1, 0, y_t],
                              [0, 0, 1, z_t],
                              [0, 0, 0, 1]])
                          
    # format df for matrix algebra                          
    trans_vecs = df_MET[['x', 'y', 'z', 'w']] 
    trans_vecs = trans_vecs.transpose()
    trans_vecs.columns = df_MET['ATOM'].tolist()      

    # multiply the methionine atom vectors by homogeneous transform matrix
    trans_vecs.CG = trans_matrix * np.matrix(trans_vecs.CG).transpose() 
    trans_vecs.SD = trans_matrix * np.matrix(trans_vecs.SD).transpose() 
    trans_vecs.CE = trans_matrix * np.matrix(trans_vecs.CE).transpose()    

    # format back for output
    trans_vecs = trans_vecs.transpose().reset_index()
    df_translated = df_MET[['CHAIN', 'RES', 'RES NAME']].reset_index(drop=True)
    df_translated = pd.concat([df_translated, trans_vecs], axis=1)
    df_translated = df_translated.rename(columns={'index': 'ATOM'})
    df_translated = df_translated.drop(['w'], axis=1)
                                       
    # append to aromatic coordinates
    df_coords = pd.concat((df_ARO, df_translated), ignore_index=True)   
    df_coords['ATOM'] = df_coords['ATOM'].str.slice(0, 1)                

    # generate xyz file
    filename = r'orca_workdir/trans_{}_{}'.format(PDB_CODE, INCR_TRANS)
    get_XYZ_file(df_coords[['ATOM', 'x', 'y', 'z']].values, filename)  # works perfectly!









