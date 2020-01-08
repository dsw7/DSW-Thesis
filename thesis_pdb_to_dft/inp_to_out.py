"""
This script creates a shell script for executing all orca operations.
"""

import os

# scope out pwd to find .xyz files
local_dirs = os.listdir('orca_workdir')

# generate input and output deck names
files_input = [i for i in local_dirs if i[-4:] == '.inp']
files_output = [i.replace('.inp', '.out') for i in files_input]

in_str = ''
for f in files_input:
    in_str += f + '\n'

array_input = """
inp_decks=({})
""".format(in_str)

out_str = ''
for f in files_output:
    out_str += f + '\n'

array_output = """
out_decks=({})
""".format(out_str)

loop = """
for((i=0; i<={}; i++));
do
    orca ${{inp_decks[i]}} > ${{out_decks[i]}}
done
""".format(len(files_input))

# merge all code
output_code = array_input + array_output + loop

# write to a shell file
with open(r'orca_workdir/inp_to_out.sh', 'w') as shellfile:
    shellfile.write(output_code)



    

