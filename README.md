# DSW-Thesis   

All code written by David Weber as part of his M.Sc. graduate work.  

---   
## Description of directories:

    ~/DSW-Thesis/thesis_Bridging_workdir/                // contains high level Met-aromatic implementation for getting bridges
    ~/DSW-Thesis/thesis_Chain_comparator/                // contains all code for assessing bridge / aromatic chain superimposition
    ~/DSW-Thesis/thesis_Met-aromatic_workdir_LOWLEVEL    // contains a lower level implementation of Met-aromatic algorithm
    ~/DSW-Thesis/thesis_Met-aromatic_workdir_HIGHLEVEL   // contains a high level Pandas/BioPython Met-aromatic algorithm implementation
    ~/DSW-Thesis/thesis_Met-aromatic_cplusplus           // contains a very basic C++ implementation of Met-aromatic algorithm
    ~/DSW-Thesis/thesis_PDB_to_DFT                       // contains workflow for translating DMS approx. relative to aromatic plane + DFT
    ~/DSW-Thesis/thesis_SerialIO                         // contains example MATLAB code for Warren Lab laser hardware serial communication
    ~/DSW-Thesis/thesis_PDB_to_DFT_part2                 // second workflow for computing Morse surfaces as a function of ||v||, Met-theta angle
    ~/DSW-Thesis/xyz_file_generator.py                   // contains a homemade .xyz file generator
    ~/DSW-Thesis/PyOpenGL_MetAromatic_Animation          // contains a PyOpenGL/pygame animation of the Met-aromatic algorithm
    ~/DSW-Thesis/legacy_Chain_comparator.tar             // legacy work for chain membership study of 2,611 structures*
    ~/DSW-Thesis/legacy_Chain_comparator_2.tar           // legacy work for chain membership study of entire PDB*
    ~/DSW-Thesis/thesis_Met-aromatic_workdir_LOWLEVEL_2  // contains a lower level implementation of Met-aromatic algorithm - REFACTORED SUMMER 2019

    *These legacy directories may be removed in 2019 as this study was done multiple times

---   
## Chain comparator sequence of events:

    (1) ::/legacy_Chain_comparator.tar           -> 2,611 proteins (Summer 2018)
    (2) ::/legacy_Chain_comparator_2.tar         -> entire PDB / IS condition assessed FOR A:MET:B (A & B = Tyr || Trp) ONLY (Sept 2018)
    (3) ::/BridgingInteractions/superimposition/ -> entire PDB / IS condition assessed FOR A:MET:B (A || B = Tyr || Trp) (Nov 2018)

