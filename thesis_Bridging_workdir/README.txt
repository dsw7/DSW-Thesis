* David *

TODO: replace with both a high level and low level implementation?

This folder contains an elegant implementation of the Met-aromatic algorithm but we go
one step further: here we compute the bridging interaction.

# --------------------------------------------------------------------------------------
Instructions:

(1) Open the script in a Python editor.
(2) Input the PDB code of interest by assigning the PDB code to PDB_CODE.
(3) Input the cutoff distance by assigning the distance (in Angstroms) to D_C.
(4) Input the cutoff angle by assigning the angle (in degrees) to A_C.
(5) Compile.
(6) Data will be outputted to the console in the following form:

    Bridges:
    Rodrigues method:
        ARO  RES  MET  BRDG
    0   PHE   38   36     3
    4   PHE  255   36     3
    8   TYR    6   36     3
    9   PHE   87  101     5
    12  PHE  104  101     5
    13  PHE  114  101     5
    19  TYR   96  101     5
    20  TYR  118  101     5
    24  PHE   87  122     2
    30  PHE  125  122     2
    Cross product method:
        ARO  RES  MET  BRDG
    0   PHE   38   36     3
    4   PHE  255   36     3
    8   TYR    6   36     3
    9   PHE   87  101     5
    12  PHE  104  101     5
    13  PHE  114  101     5
    19  TYR   96  101     5
    20  TYR  118  101     5
    24  PHE   87  122     3
    30  PHE  125  122     3
    34  TYR  118  122     3
    
(7) BRDG returns how many aromatic residues meet Met-aromatic criteria relative to a methionine residue MET.