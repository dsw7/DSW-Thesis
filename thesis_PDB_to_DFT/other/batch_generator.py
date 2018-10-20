"""
My first batch/shell writer
Worked last time I used it
May not necessarily include in thesis but it's good for getting Windows .bat generator code
I used mostly MacOS for DFT pretty much because I only have ORCA installed on school Windows platform
"""

import os
from platform import system

# scope out local environment and find .inp decks
workdir = os.getcwd()
listdir = os.listdir(workdir)
dirs_input = []
for dirs in listdir:
    _, file_extension = os.path.splitext(dirs)
    if file_extension == '.inp':
        dirs_input.append(dirs)



if system() == 'Windows':            
    # create actual file        
    file_name = 'ORCA_batchprocessor'
    file = open(file_name + '.bat', 'w')   
    
    file.write("@echo off\n")    
    file.write("setlocal enabledelayedexpansion\n")
    cd = 'cd {}\n'.format(workdir)
    file.write(cd)
    file.write('\n')
    
    # generate input arguments
    for index, filename in enumerate(dirs_input):
        line='set in[{a}]={b}\n'.format(a=index, b=filename)
        file.write(line)    
    file.write('\n')
    
    # generate output arguments
    for index, filename in enumerate(dirs_input):
        line = 'set out[{a}]={b}.out\n'.format(a=index, b=filename[:-4])
        file.write(line)
    file.write('\n')
    
    # write the actual loop
    loop = 'for /l %%n in (0,1,{}) do (\n'.format(len(dirs_input))
    file.write(loop)   
    file.write("   IF EXIST !in[%%n]! (\n")
    file.write("   echo Processing !in[%%n]!\n")
    file.write("   orca !in[%%n]! > !out[%%n]!\n")
    file.write("   ) ELSE (\n")
    file.write("   echo ORCA file does not exist.\n")
    file.write("   )\n")
    file.write(")\n")
    file.write("PAUSE\n")
    
    file.close()

else:    
    # create actual file        
    file_name = 'ORCA_batchprocessor'
    file = open(file_name + '.sh', 'w')   
    
    file.write("#!/bin/bash\n")    
    file.write("#type ./ORCA_batchprocessor.sh into terminal to run script\n")
    file.write("#type in the event of permission error: sudo chmod 755 'ORCA_batchprocessor.sh'\n")
    cd = 'cd {}\n'.format(workdir)
    file.write(cd)
    file.write('\n')
    
    # generate input shell array
    a_line = 'inputs=('
    s = str(dirs_input).replace(',', '')
    s = s.replace("'", '')
    s = s[1:-1]
    line = a_line + s + ')\n'
    file.write(line)    
    
    # generate output shell array
    a_line = 'outputs=('
    s = s.replace(".inp", '.out')
    line = a_line + s + ')\n'
    file.write(line)    

    file.write('for ((i=0;i<${#inputs[@]};++i)); do\n')
    file.write('    SECONDS=0\n')
    file.write('    orca ${inputs[i]} > ${outputs[i]};\n')
    file.write('    duration=$SECONDS\n')
    file.write('    echo $duration\n')
    file.write('done\n')
    file.write('echo End\n')
    
    file.close()
    
    
