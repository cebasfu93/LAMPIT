InputCore = au144SR60.pdb
InputLigand = ligand.mol2
OutputFile = Au_144-SN_60.pdb
ResidueName = LFU
AnchorIndexes = 0,6,12

VMD := /software/vmd-1.9.1/bin/vmd191

$(OutputFile) : $(InputCore) $(InputLigand)
	python3.6 NP_builder.py -l $(InputLigand) -c $(InputCore) -o $(OutputFile) -r $(ResidueName) -a $(AnchorIndexes)
	$(VMD) $(OutputFile) $(InputLigand)
	#cp $(OutputFile) ../

clean:
	rm -rf $(OutputFile)
