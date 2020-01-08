"""
Use this script to get isolate a single Met-aromatic interaction. Can be used
for single interaction DFT calculations.
"""

PDB_CODE = '1rcy'  # input the PDB code of interest
AROMATIC_RESIDUE = 122
METHIONINE_RESIDUE = 18


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# dependencies

import sys
import Bio.PDB as bp
import os
import pandas as pd

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
    
    
# merge both dfs / format atom
df_coords = pd.concat((df_ARO, df_MET))
df_coords['ATOM'] = df_coords['ATOM'].str.slice(0, 1)

# export data
get_XYZ_file(df_coords[['ATOM', 'x', 'y', 'z']].values, 'myfile')  # works perfectly!










