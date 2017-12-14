#!/bin/bash

#OPTIONS FOR THE USER
SYS_NAME="NP1-1"  #Prefix used for the generation of files
LIG_MOL2="Mol-ia1_m1-c2" #mol2 file of a linear ligand with the right charges
CORE_PDB="au144SR60"  #path to the pdb of the NP's core including the first carbon
ANCHOR_NDX="0,6,12,18,24,30,34,40,44,48,54,58,62,68,72"  #Indexes of the atoms in LIG_MOL2 that will be aligned to the COM-C1 vector
OLD_NAME="F00"  #Name of the ligand in LIG_MOL2
NEW_NAME="LF1"  #New 3-letter residue name for the coating

F_LEAP1="LeapLig"  #Name of the first tleap input
F_LEAP2="LeapSys"  #Name of the second tleap input
DEPENDS="/DATA/SoftwareSFU/IN-HOUSE/LAMPIT/DEPENDENCIES"  #Path of the folder with LAMPIT's dependencies
ANCHOR_NAME="C1"
ANCHOR_H="H1,H2"

#Creates working directory and copies important files
mkdir ${SYS_NAME}
cp ${LIG_MOL2}.mol2 ${CORE_PDB}.pdb ${SYS_NAME}
cp LAMPIT.sh ${SYS_NAME}/LAMPIT-run.sh

#CHANGES RESIDUE NAME AND MAKES NEW MOL3 FILE
sed  s"/${OLD_NAME}/${NEW_NAME}/" ${SYS_NAME}/${LIG_MOL2}.mol2 > ${SYS_NAME}/${NEW_NAME}.mol2

#STONES THE SPECIFIED ATOMS OF THE LIGAND TO THE SULPHURS OF THE NP
python3.6 ${DEPENDS}/NP_builder.py -l ${SYS_NAME}/${NEW_NAME}.mol2 -c ${SYS_NAME}/${CORE_PDB}.pdb -o ${SYS_NAME}/${SYS_NAME}_stoned.pdb -r ${NEW_NAME} -s ${ANCHOR_NDX}

#WRITED FILE WITH THE FF PARAMETERS FOR THE LIGAND
parmchk2 -i ${SYS_NAME}/${NEW_NAME}.mol2 -f mol2 -o ${SYS_NAME}/${NEW_NAME}.frcmod -a y

#Writes first input file for tleap and writes the .lib file for the ligand
echo "source leaprc.gaff" > ${SYS_NAME}/${F_LEAP1}.in
echo "loadamberparams ${NEW_NAME}.frcmod \n" >> ${SYS_NAME}/${F_LEAP1}.in
echo "${NEW_NAME} = loadmol3 ${NEW_NAME}.mol2" >> ${SYS_NAME}/${F_LEAP1}.in
echo "check ${NEW_NAME}" >> ${SYS_NAME}/${F_LEAP1}.in
echo "saveoff ${NEW_NAME} ${NEW_NAME}.lib \n" >> ${SYS_NAME}/${F_LEAP1}.in
echo "#saveamberparm ${NEW_NAME} ${NEW_NAME}.prmtop ${NEW_NAME}.inpcrd \n" >> ${SYS_NAME}/${F_LEAP1}.in
echo "quit" >> ${SYS_NAME}/${F_LEAP1}.in

tleap -sf ${SYS_NAME}/$F_LEAP1.in > ${SYS_NAME}/${F_LEAP1}.log

#Writes second input file for tleap the .prmtop and .inpcrd for the complete system
echo "source leaprc.gaff \n" > ${SYS_NAME}/${F_LEAP2}.in
echo "loadoff ${SYS_NAME}/${NEW_NAME}.lib" >> ${SYS_NAME}/${F_LEAP2}.in
echo "loadamberparams ${SYS_NAME}/${NEW_NAME}.frcmod" >> ${SYS_NAME}/${F_LEAP2}.in
echo "loadamberparams ${DEPENDS}/AU.frcmod" >> ${SYS_NAME}/${F_LEAP2}.in
echo "loadamberparams ${DEPENDS}/ST.frcmod \n" >> ${SYS_NAME}/${F_LEAP2}.in
echo "${NEW_NAME} = loadmol3 ${SYS_NAME}/${NEW_NAME}.mol2" >> ${SYS_NAME}/${F_LEAP2}.in
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

#Modifies topology file to include bonds and angles involving staple atoms
python3.6 ${DEPENDS}/staples_topology.py -p ${SYS_NAME}/${SYS_NAME}.top -x ${SYS_NAME}/${SYS_NAME}.gro -a ${ANCHOR_NAME} -y ${ANCHOR_H} -f ${SYS_NAME}
