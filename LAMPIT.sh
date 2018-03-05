#!/bin/bash

#OPTIONS FOR THE USER
SYS_NAME="test"  #Prefix used for the generation of files
CORE_PDB="au144SR60"  #path to the pdb of the NP's core including the first carbon

MOL2_LIG1="LF2" #mol2 file of a linear ligand 1 with the right charges
OLD_NAME1="F00"  #Name of the ligand 1 in MOL2_LIG1
NEW_NAME1="LF2"  #New 3-letter residue name for ligand 1

MOL2_LIG2="LF3" #mol2 file of a linear ligand 2 with the right charges
OLD_NAME2="F00"  #Name of the ligand 2 in MOL2_LIG1
NEW_NAME2="LF3"  #New 3-letter residue name for ligand 2

FRAC_LIG1="1.0"

CORENAME="Au"  #Name core atom in CORE_PDB
STAPLENAME="S" #N
COREANCHOR="C"
MORPHOLOGY="random"
RSEED="666"

F_LEAP1="LeapLig"  #Name of the first tleap input
F_LEAP2="LeapSys"  #Name of the second tleap input
F_NPBUILDER="NP_builder"
F_STAPLES="staples"
DEPENDS="/DATA/SoftwareSFU/IN-HOUSE/LAMPIT/DEPENDENCIES/"  #Path of the folder with LAMPIT's dependencies

########################################################################################################################################################################################
#Creates working directory and copies important files
mkdir ${SYS_NAME}
cp ${MOL2_LIG1}.mol2 ${SYS_NAME}/${NEW_NAME1}.mol2
cp ${MOL2_LIG2}.mol2 ${SYS_NAME}/${NEW_NAME2}.mol2
cp ${CORE_PDB}.pdb ${SYS_NAME}/
cp LAMPIT.sh ${SYS_NAME}/LAMPIT-run.sh

#CHANGES RESIDUE NAME AND MAKES NEW MOL3 FILE
sed -i s"/${OLD_NAME1}/${NEW_NAME1}/" ${SYS_NAME}/${NEW_NAME1}.mol2
sed -i s"/${OLD_NAME1}/${NEW_NAME2}/" ${SYS_NAME}/${NEW_NAME2}.mol2

echo "output \t ${SYS_NAME}/${SYS_NAME}_stoned.pdb" > ${SYS_NAME}/${F_NPBUILDER}.in
echo "core \t ${SYS_NAME}/${CORE_PDB}.pdb" >> ${SYS_NAME}/${F_NPBUILDER}.in
echo "corename \t ${CORENAME}" >> ${SYS_NAME}/${F_NPBUILDER}.in
echo "staplename \t ${STAPLENAME}" >> ${SYS_NAME}/${F_NPBUILDER}.in
echo "coreanchor \t ${COREANCHOR}" >> ${SYS_NAME}/${F_NPBUILDER}.in
echo "ligand1 \t ${SYS_NAME}/${NEW_NAME1}.mol2" >> ${SYS_NAME}/${F_NPBUILDER}.in
echo "ligand2 \t ${SYS_NAME}/${NEW_NAME2}.mol2" >> ${SYS_NAME}/${F_NPBUILDER}.in
echo "morphology \t ${MORPHOLOGY}" >> ${SYS_NAME}/${F_NPBUILDER}.in
echo "frac_lig1 \t ${FRAC_LIG1}" >> ${SYS_NAME}/${F_NPBUILDER}.in
echo "rseed \t ${RSEED}" >> ${SYS_NAME}/${F_NPBUILDER}.in

#STONES THE SPECIFIED ATOMS OF THE LIGAND TO THE SULPHURS OF THE NP
python3.6 ${DEPENDS}/NP_builder.py ${SYS_NAME}/${F_NPBUILDER}.in


#WRITED FILE WITH THE FF PARAMETERS FOR THE LIGAND
parmchk2 -i ${SYS_NAME}/${NEW_NAME1}.mol2 -f mol2 -o ${SYS_NAME}/${NEW_NAME1}.frcmod -a y
parmchk2 -i ${SYS_NAME}/${NEW_NAME2}.mol2 -f mol2 -o ${SYS_NAME}/${NEW_NAME2}.frcmod -a y

#Writes first input file for tleap and writes the .lib file for the ligand
echo "source leaprc.gaff" > ${SYS_NAME}/${F_LEAP1}.in
echo "loadamberparams ${SYS_NAME}/${NEW_NAME1}.frcmod \n" >> ${SYS_NAME}/${F_LEAP1}.in
echo "loadamberparams ${SYS_NAME}/${NEW_NAME2}.frcmod \n" >> ${SYS_NAME}/${F_LEAP1}.in
echo "${NEW_NAME1} = loadmol3 ${SYS_NAME}/${NEW_NAME1}.mol2" >> ${SYS_NAME}/${F_LEAP1}.in
echo "${NEW_NAME2} = loadmol3 ${SYS_NAME}/${NEW_NAME2}.mol2" >> ${SYS_NAME}/${F_LEAP1}.in
echo "check ${NEW_NAME1}" >> ${SYS_NAME}/${F_LEAP1}.in
echo "check ${NEW_NAME2}" >> ${SYS_NAME}/${F_LEAP2}.in
echo "saveoff ${NEW_NAME1} ${SYS_NAME}/${NEW_NAME1}.lib \n" >> ${SYS_NAME}/${F_LEAP1}.in
echo "saveoff ${NEW_NAME2} ${SYS_NAME}/${NEW_NAME2}.lib \n" >> ${SYS_NAME}/${F_LEAP1}.in
echo "quit" >> ${SYS_NAME}/${F_LEAP1}.in

tleap -sf ${SYS_NAME}/$F_LEAP1.in > ${SYS_NAME}/${F_LEAP1}.log

#Writes second input file for tleap the .prmtop and .inpcrd for the complete system
echo "source leaprc.gaff \n" > ${SYS_NAME}/${F_LEAP2}.in
echo "loadoff ${SYS_NAME}/${NEW_NAME1}.lib" >> ${SYS_NAME}/${F_LEAP2}.in
echo "loadoff ${SYS_NAME}/${NEW_NAME2}.lib" >> ${SYS_NAME}/${F_LEAP2}.in
echo "loadamberparams ${SYS_NAME}/${NEW_NAME1}.frcmod" >> ${SYS_NAME}/${F_LEAP2}.in
echo "loadamberparams ${SYS_NAME}/${NEW_NAME2}.frcmod" >> ${SYS_NAME}/${F_LEAP2}.in
echo "loadamberparams ${DEPENDS}/AU.frcmod" >> ${SYS_NAME}/${F_LEAP2}.in
echo "loadamberparams ${DEPENDS}/ST.frcmod \n" >> ${SYS_NAME}/${F_LEAP2}.in
echo "${NEW_NAME1} = loadmol3 ${SYS_NAME}/${NEW_NAME1}.mol2" >> ${SYS_NAME}/${F_LEAP2}.in
echo "${NEW_NAME2} = loadmol3 ${SYS_NAME}/${NEW_NAME2}.mol2" >> ${SYS_NAME}/${F_LEAP2}.in
echo "AU = loadmol3 ${DEPENDS}/AU.mol2" >> ${SYS_NAME}/${F_LEAP2}.in
echo "ST = loadmol3 ${DEPENDS}/ST.mol2 \n" >> ${SYS_NAME}/${F_LEAP2}.in
echo "${SYS_NAME} = loadpdb ${SYS_NAME}/${SYS_NAME}_stoned.pdb"  >> ${SYS_NAME}/${F_LEAP2}.in
echo "savepdb ${SYS_NAME} ${SYS_NAME}/${SYS_NAME}.pdb \n" >> ${SYS_NAME}/${F_LEAP2}.in

echo "saveamberparm ${SYS_NAME} ${SYS_NAME}/${SYS_NAME}.prmtop ${SYS_NAME}/${SYS_NAME}.inpcrd \n" >> ${SYS_NAME}/${F_LEAP2}.in
echo "quit" >> ${SYS_NAME}/${F_LEAP2}.in

tleap -sf ${SYS_NAME}/${F_LEAP2}.in > ${SYS_NAME}/${F_LEAP2}.log

#Converts amber geometry and topology to gromacs format
python ${DEPENDS}/acpype.py -p ${SYS_NAME}/${SYS_NAME}.prmtop -x ${SYS_NAME}/${SYS_NAME}.inpcrd -b ${SYS_NAME} -c user -r > acpype.log
mv acpype.log ${SYS_NAME}
mv ${SYS_NAME}_GMX.gro ${SYS_NAME}/${SYS_NAME}.gro
mv ${SYS_NAME}_GMX.top ${SYS_NAME}/${SYS_NAME}.top

#Remove files that ABSOLUTELY unnecessary
rm -rf em.mdp md.mdp leap.log ANTECHAMBER.FRCMOD


echo "topfile \t ${SYS_NAME}/${SYS_NAME}.top" > ${SYS_NAME}/${F_STAPLES}.in
echo "grofile \t ${SYS_NAME}/${SYS_NAME}.gro" >> ${SYS_NAME}/${F_STAPLES}.in
echo "resname1 ${NEW_NAME1}\t " >> ${SYS_NAME}/${F_STAPLES}.in
echo "resname2 ${NEW_NAME2}\t " >> ${SYS_NAME}/${F_STAPLES}.in
echo "ligand1 \t ${SYS_NAME}/${NEW_NAME1}.mol2" >> ${SYS_NAME}/${F_STAPLES}.in
echo "ligand2 \t ${SYS_NAME}/${NEW_NAME2}.mol2" >> ${SYS_NAME}/${F_STAPLES}.in
echo "workdir \t ${SYS_NAME}" >> ${SYS_NAME}/${F_STAPLES}.in


#Modifies topology file to include bonds and angles involving staple atoms
python3.6 ${DEPENDS}/staples.py ${SYS_NAME}/${F_STAPLES}.in
