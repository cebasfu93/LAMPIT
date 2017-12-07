#!/bin/bash

#OPTIONS FOR THE USER
LIG_MOL2="Mol-ia1_m1-c2.mol2"
CORE_PDB="au144SR60.pdb"
ANCHOR_NDX ="0,6,12"
OLD_NAME="F00"

NEW_NAME="LF1"
STONED_PDB="stoned-test.pdb"
F_LEAP1="LeapLig.in"
F_LEAP2="LeapSys.in"
F_GMX="sys"
SYS_PDB="./"

#CHANGES RESIDUE NAME AND MAKES NEW MOL3 FILE
sed  s'/$OLD_NAME/$NEW_NAME/' $LIG_MOL2 > $NEW_NAME.mol2

#STONES THE SPECIFIED ATOMS OF THE LIGAND TO THE SULPHURS OF THE NP
python3.6 NP_builder -l $NEW_NAME.mol2 -c $CORE_PDB -o $STONED_PDB -r $NEW_NAME -a $ANCHOR_NDX

#WRITED FILE WITH THE FF PARAMETERS FOR THE LIGAND
parmchk2 -i $NEW_NAME.mol2 -f mol2 -o $NEW_NAME.frcmod -a y

#
echo "source leaprc.gaff \n" > $F_LEAP1
echo "loadamberparms $NEW_NAME.frcmod \n" >> $F_LEAP1
echo "$NEW_NAME = loadmol2 $NEW_NAME.mol2 \n" >> $F_LEAP1
echo "check $NEW_NAME \n" >> $F_LEAP1
echo "saveoff $NEW_NAME $NEW_NAME.lib \n" >> $F_LEAP1
echo "#saveamberparm $NEW_NAME $NEW_NAME.prmtop $NEW_NAME.inpcrd \n" >> $F_LEAP1
echo "quit" >> $F_LEAP1

#
tleap -sf $F_LEAP1

#
echo "source leaprc.gaff \n" > $F_LEAP2

echo "loadoff $NEW_NAME.lib \n" >> $F_LEAP2

echo "loadamberparams $NEW_NAME.frcmod \n" >> $F_LEAP2
echo "loadamberparams AU.frcmod \n" >> $F_LEAP2
echo "loadamberparams ST.frcmod \n" >> $F_LEAP2

echo "$NEW_NAME = loadmol3 $NEW_NAME.mol2 \n" >> $F_LEAP2
echo "AU = loadmol3 AU.mol3 \n" >> $F_LEAP2
echo "ST = loadmol3 ST.mol3 \n" >> $F_LEAP2
echo "system = loadpdb $SYS_PDB \n" >> $F_LEAP2

echo "saveamberparm system $F_GMX.prmtop $F_GMX.inpcrd \n" >> $F_LEAP2
echo "quit"

#
tleap -sf $F_LEAP2

#
acpype -p $F_GMX.prmtop -x $F_GMX.inpcrd -r
