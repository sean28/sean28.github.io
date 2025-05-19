### File format converter
<div style="text-align: justify">The conversion between various chemical formats is a headache question. Different chemical formats have different uses.For example, In the field of molecular dynamics simulation, the format is used as follows:</div>


|Name|Usage|
|:---|:---|
|Coordinate file.     |The coordinate file records the three-dimensional coordinates of all atoms in the simulation system. Format include: amber (prmtop, pram7) gromacs (top) |
|Topology file.   |Topological files record the connection relations of atoms and molecular mechanical parameters. Format include: amber (inpcrd) gromacs (gro)|
|Trajectory data    |Record the coordinate file of each frame atom. Format include: gromacs (trr, xtc) amber (mdcrd) namd(dcd)|
|Reference structure|Atomic coordinate files used as reference structures in some software and sometimes as topology files. The Reference files needs to be converted into topology file and coordinate file for calculation. Format include: pdb, gro, mol2, sdf|


<div style="text-align: justify"> Different software supports different formats. The same format converted by the same software may be different.When selecting the reference and trajectory files, they must be the reference structure files and trajectory files saved after the same software aligned, otherwise the location information error will occur! Therefore, it is very important to skillfully use a format chemical format conversion tool.</div>
<div style="text-align: justify"> <br> </div>
<div style="text-align: justify"> Open Babel is a chemical toolbox designed to speak the many languages of chemical data. It’s an open, collaborative project allowing anyone to search, convert, analyze, or store data from molecular modeling, chemistry, solid-state materials, biochemistry, or related areas.However, it is difficult for beginners to get started with this software. <a href="http://www.cheminfo.org/Chemistry/Cheminformatics/FormatConverter/index.html">Here</a> is an online openbabel tool for easy to use. </div>
<div style="text-align: justify"> <br> </div>
<div style="text-align: justify"> In addition, a script for converting smiles numbers into SDF in batch is also provided  <a href="https://drive.google.com/file/d/1k_pAnCFhXI2teUd5u40vghfXk2E51a3a/view?usp=sharing">here</a>. This script needs an excel table, and the corresponding smiles number is filled in.</div>
<div style="text-align: justify"> <br> </div>
<div style="text-align: justify"> The table format is as follows: </div>

|Compound name|SMILES|
|:---|:---|
|1 |CCCCOc1ccccc1|
|2 |CCCCOc1ccccc1|
|3 |CCCCOc1ccccc1|
|..|..           |

The scripts usage is as follows: 

```
python python smiles_to_2dsdf.py xx.xlsx
```