LAMPIT: Look At Me Putting Interesting Names

Script in bash/python for the preparation of a NP coated by a given ligand.
From xyz of NP core and mol2 of the ligand, the results are the parameter and coordinates of the coated NP.


Core (.pdb)  Ligand (.mol2)  ----> Coated NP ----> Coordinates (.gro)  Topology (.top)

Remarks:
1. Usage: Download LAMPIT.sh and DEPENDENCIES/, then modify the variables in LAMPIT.sh and run.
2. To use it with only one ligand, the option FRAC_LIG1 must be set to 1.0.
4. When using 2 ligands, random and janus distributions are supported.
5. Ligand mol2 file must have the keywords @<TRIPOS>ATOM, @<TRIPOS>BOND, and @<TRIPOS>RESIDUECONNECT.
6. Immediately beneath @<TRIPOS>RESIDUECONNECT the first written thing must be the number of the atom used as anchor.
7. The output name for the core atoms will always be "AU"
8. The output name for the stapling atom will always be "ST"
9. The options in NP_builder.in and staples.in are case sensitive
