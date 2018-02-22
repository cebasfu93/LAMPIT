LAMPIT: Look At Me Putting Interesting Names

Script in bash/python for the preparation of a NP coated by a given ligand.
From xyz of NP core and mol2 of the ligand, the results are the parameter and coordinates of the coated NP.


Core (.xyz)	  		     Coordinates (.gro)
	       ----> Coated NP ---->
Ligand (.mol2)			     Topology (.top)

Remarks:
1. Ligand mol2 file must have the keywords @<TRIPOS>ATOM, @<TRIPOS>BOND, and @<TRIPOS>RESIDUECONNECT.
2. Immediately beneath @<TRIPOS>RESIDUECONNECT the first written thing must be the number of the atom used as anchor.
3. The output name for the core will always be "AU"
4. The output name for the stapling atom will always be "ST"
