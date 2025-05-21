# Script for calculating aMD parameters
<div style="text-align: justify"> The time scale of traditional molecular dynamics simulation is usually hundreds of nanoseconds, and it is difficult to capture biological processes that can be observed at microsecond or even millisecond scales. Accelerated molecular dynamics (aMD) reduces the height of local energy barrier by modifying potential energy, so as to accelerate sampling.</div>
<div style="text-align: justify"> <br> </div>
<div style="text-align: justify"> When running aMD, the following four parameters need to be modified to modify the potential energy.</div>


|Parameters            |        Description          | Formula calculation         |
|          :---        |           :---              |           :---              |
|EthreshP. | Average total potential energy threshold.                     |E(tot)= EPtot (kcal/mol) + (0.16kcal/mol/atom * whole system atoms)|
|alphaP.   | Inverse strength boost factor for the total potential energy. |Alpha(tot)= (0.16 kcal/mol/atom * whole system atoms) |
|EthreshD. | Average dihedral energy threshold.                            |E(dih)= DIHED (kcal/mol) + (4 kcal/mol/residues * solute residues)|
|alphaD.   | Inverse strength boost factor for the dihedral energy.        |Alpha(dih)= 0.2 * (4 kcal/mol/residues * solute residues)|


Note: EPtot: An average total potential energy. DIHED: An average dihedral energy.

<div style="text-align: justify">A shortcut calculation script is provided here. After entering EPtot, DIHED, solute residues, and whole system atoms, the values of EtheshP, alphaP, EtheshD, and alphaD will be calculated automatically </div>


```
#!/usr/bin/python
# -*- coding:utf-8 -*-
# Created By Sean -- MUST ;
# Caculate  EthreshP alphaP EthreshD alphaD value of Amber aMD 

# input
EPtot = input("Please enter a value for EPtot(kcal/mol)：")
DIHED = input("Please enter a value for DIHED(kcal/mol)：")
atom_num = input("Please enter atoms：")
resi_num = input("Please enter solute residues：")

EPtot = float(EPtot)
DIHED = float(DIHED)

# caculate 
EthreshP =  EPtot + (0.16 * atom_num);
alphaP = 0.16 * atom_num;
EthreshD = DIHED + (4 * resi_num);
alphaD = 0.2 * (4 * resi_num);

EthreshP = round(EthreshP, 2)
alphaP = round(alphaP, 2)
EthreshD = round(EthreshD, 2)
alphaD = round(alphaD, 2)

# output 
print ("EthreshP = ", EthreshP);
print ("alphaP = ", alphaP);
print ("EthreshD = ", EthreshD);
print ("alphaD = ", alphaD);

```
