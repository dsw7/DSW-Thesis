# DSW Automated Pipeline for Performing ORCA Density Functional Theory Calculations

## Description  
This is a super interesting directory to both chemists and computer scientists alike. The directory
gets an aromatic interaction of choosing from the PDB. A homogeneous transformation matrix
is then used to increment the position of a dimethyl sulfide molecule (approximation for the
methionine CE-SD-CG moiety) along a vector normal to the aromatic plane. Here we can force
the DMS molecule closer or further away from the aromatic ring. These translations are then
exported using a homemade `.xyz` file generator for rendering hydrogens in Avogadro. Then a Python
script writes homemade `.inp` files for ORCA. The most interesting part is that there
is a script that writes code, "code that writes code", in a process analogous to a computer virus.
The Python script writes a shell script for the execution of DFT calculations using the premade
input decks. This shell script then executes all DFT calculations, creating `.out` decks in the process.
Lastly, a script is used to read the `.out` decks and plot the HOMO and LUMO information for each
translation as a function of the translation.


---  
## Instructions:

### Step 1

Open `get_translated.py` and select the methionine and aromatic residue position numbers
corresponding to a Met-aromatic interaction by manually editing AROMATIC_RESIDUE and METHIONINE_RESIDUE.
Run the script. This will place all the `.xyz` files into:

    ~/thesis_pdb_to_dft/orca_workdir

### Step 2

Open the `.xyz` files in Avogadro and manually add hydrogens. Save the files. Yes... I know this is frustrating.

### Step 3

Next we generate input decks using `xyz_to_inp.py`. Simply run the script.

### Step 4

Then we generate the shell script for running ORCA calculations using `inp_to_out.py`. 
Simply run the script. The shell script `inp_to_out.sh` will be located in:

    ~/thesis_pdb_to_dft/orca_workdir

### Step 5

Run the shell script. To do so, open up cmd or terminal, depending on OS, then `cd` into `orca_workdir`.
For me, it is:  
    
    $ cd ~/thesis_pdb_to_dft/orca_workdir

Then type:   

    $ sh inp_to_out.sh  

And ENTER. The DFT calculation will now begin. I used the default DFT parameters which render the DFT calculations very quick. Be advised that the amount of time taken to complete the DFT will increase significantly if my code is modified with other functionals. Increase your fan speeds if you do change my functionals as your computer may get super hot.

### Step 6

If all goes well, `~/thesis_pdb_to_dft/orca_workdir` will be loaded with many `.out` files. Use cleanup.py to clean out
other random junk files, unless you need them. This will help keep `~/orca_workdir` clean. DFT calculations are notorious
for generating a ton of random files. Note that you may still need to go through and clean out the folder manually.

### Step 7

Use `get_energies.py` to generate a plot for the translation series. The x-axis may need to be flipped.

### Step 8

I highly recommend completely clearing `~/thesis_pdb_to_dft/orca_workdir` as soon as the experiment is done.
This folder will get filled with lots of crap and fast.

---   
## Other notes:  

    ~/thesis_pdb_to_dft/get_raw.py // can be used to fetch a single, untranslated, and non-geometry optimized "met-aromatic" pair
