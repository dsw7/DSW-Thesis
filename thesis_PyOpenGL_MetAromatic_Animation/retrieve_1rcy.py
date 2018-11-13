# function for importing raw .pdb coordinates

IDX_ATOM = 0
IDX_CHAIN = 4
CHAIN = 'A'

def import_coords():
    # retrieve working coordinates from PDB file 1rcy
    with open('1rcy.pdb', 'r') as f:
        data_incoming = f.readlines()

    # chunk lines 
    data = [line.split() for line in data_incoming]

    # stop at end of first model - fix for bug in v1.3
    model_first = []
    for line in data:
        if line[IDX_ATOM] != 'ENDMDL':
            model_first.append(line)
        else:
            break
        
    # get only ATOM records
    model_first = [line for line in model_first if line[IDX_ATOM] == 'ATOM']

    # get only specific chains
    return [line for line in model_first if line[IDX_CHAIN] == CHAIN]


