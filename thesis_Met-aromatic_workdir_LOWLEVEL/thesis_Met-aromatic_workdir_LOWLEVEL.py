# example usage of met-aromatic low level routine

from sys                    import path as PATH
from os                     import getcwd, path
PATH.append(path.join(getcwd(), 'lib'))

from PDB_filegetter         import PDBFile
from ma_lowlevel            import met_aromatic

CODE = '4n07'     # input some PDB code here
CHAIN = 'A'       # input a chain of interest here
CUTOFF = 6.0      # input some cutoff vector v norm
ANGLE = 109.5     # input some cutoff vector a / vector v and vector g / vector v angle
MODEL = 'cp'      # 'cp' or 'rm' -> Cross Product or Rodrigues' Method lone pair interpolation methods

file_pdb = PDBFile(CODE)
path_to_file = file_pdb.fetch_from_PDB()

data_retrieved = met_aromatic(CHAIN=CHAIN, CUTOFF=CUTOFF, ANGLE=ANGLE, 
                              MODEL=MODEL, filepath=path_to_file)
        
for item in data_retrieved:  # do something with the data - here I print to console
    print(item)

file_pdb.clear()
