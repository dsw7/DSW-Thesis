# example usage of met-aromatic low level routine

from sys                    import path as PATH
from os                     import getcwd, path
PATH.append(path.join(getcwd(), 'lib'))

from PDB_filegetter         import PDBFile
from ma_lowlevel            import met_aromatic

CODE = '4n07'
CHAIN = 'A'
CUTOFF = 6.0
ANGLE = 109.5
MODEL = 'cp'

file_pdb = PDBFile(CODE)
path_to_file = file_pdb.fetch_from_PDB()

data_retrieved = met_aromatic(CHAIN=CHAIN, CUTOFF=CUTOFF, ANGLE=ANGLE, 
                              MODEL=MODEL, filepath=path_to_file)
        
for item in data_retrieved:
    print(item)

file_pdb.clear()