"""
This script converts Avogadro .xyz files into .inp input decks.
"""

import os

# scope out pwd to find .xyz files
local_dirs = os.listdir('orca_workdir')
local_dirs = [i for i in local_dirs if i[-4:] == '.xyz']  # or use os.splitext()

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
    with open(os.path.join('orca_workdir', file)) as f_INPUT:
        lines = [item for item in f_INPUT if item[0] in ['C', 'S', 'H']]
        
    input_filename = file.split('.')
    output_filename = '{}.{}.inp'.format(*input_filename)
    
    with open(os.path.join('orca_workdir', output_filename), 'w') as f_OUTPUT:
        for line in HEADER + lines + FOOTER:
            f_OUTPUT.write(line)
        
    
    

