"""
This script removes .gbw and .prop garbage.
This script is not guaranteed to remove all useless garbage.
"""

import os
cwd = 'orca_workdir'

gbw_dirs = os.listdir(cwd)
gbw_dirs = [i for i in gbw_dirs if i[-4:] == '.gbw']

for item in gbw_dirs:
    os.remove(os.path.join(cwd, item))
    
prop_dirs = os.listdir(cwd)
prop_dirs = [i for i in prop_dirs if i[-5:] == '.prop']

for item in prop_dirs:
    os.remove(os.path.join(cwd, item))
    
proptxt_dirs = os.listdir(cwd)
proptxt_dirs = [i for i in proptxt_dirs if 'property.txt' in i]

for item in proptxt_dirs:
    os.remove(os.path.join(cwd, item))