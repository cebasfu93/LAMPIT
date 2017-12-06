InputCore = Au144_MD6341surface1.xyz

InputLigand = ligand.mol2
OutputFile = Au144L1-60.pdb
ResidueName = LFU

#InputLigand = MET.mol2
#OutputFile = Au144_MET.pdb
#ResidueName = MET

VMD := /software/vmd-1.9.1/bin/vmd191

$(OutputFile) : $(InputCore) $(InputLigand)
	python3.6 NP_builder.py -l $(InputLigand) -c $(InputCore) -o $(OutputFile) -r $(ResidueName)
	$(VMD) $(OutputFile) $(InputLigand)


clean:
	rm -rf $(OutputFile)
