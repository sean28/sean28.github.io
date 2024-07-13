# Welcome to team pages of Sean drug design and discovery 

## Introduction
---

![This is an image](/cadd.png)

<div style="text-align: justify"> The content of this website mainly involves computer-aided drug design, molecular simulation, computational chemistry, computational biology and bioinformatics. This site is currently used as a learning record, and the source of reprinted articles is indicated. If there is infringement, please contact to delete. </div>

## Table of contents
---

- [1. Molecular dynamics simulation](#1-molecular-dynamics-simulation)
  * [1.1 GROMACS](#11-gromacs)
    + [An auto protonation-pdb2gmx run script](#an-auto-protonation-pdb2gmx-run-script)
    + [Auto gromacs-result analysis program](#auto-gromacs-result-analysis-program)
  * [1.2 AMBER](#12-amber)
    + [Online amber tool](#online-amber-tool)
    + [Script for calculating aMD parameters](#script-for-calculating-amd-parameters)
- [2. Perl](#2-perl)
  * [2.1 Extract selected residues form protein pdb](#21-extract-selected-residues-form-protein-pdb)
  * [2.2 Extract ligand form complex protein-ligand pdb](#22-extract-ligand-form-complex-protein-ligand-pdb)
- [3. Python](#3-python)
  * [3.1 Python basics](#31-python-basics)
    + [Python reads user input files](#python-reads-user-input-files)
  * [3.2 Python processing office](#32-python-processing-office)
    + [CSV](#csv)
    + [EXCEl](#excel)
    + [PPT](#ppt)
  * [3.3 Python processing image](#33-python-processing-image)
    + [Vertical merge png](#vertical-merge-png)
    + [Horizontal merge png](#horizontal-merge-png)
- [4. R](#4-r)
  * [4.1 Parallel kmeans scripts](#41-parallel-kmeans-scripts)
  * [4.2 ggplot2](#42-ggplot2)
  * [4.3 Binder](#43-binder)
  * [4.4 TRAPP Multiple Comparison Script](#44-trapp-multiple-comparison-script)
- [5. Pymol](#5-pymol)
- [6. Free Energy](#6-free-energy)
- [7. Others](#7-others)
  * [7.1 Automatically building scientific research environment](#71-automatically-building-scientific-research-environment)
  * [7.2 Online calculator](#72-online-calculator)
  * [7.3 File format converter](#73-file-format-converter)
  * [7.4 Online chemical structure drawing and editing tool](#74-online-chemical-structure-drawing-and-editing-tool)
- [Data Availability Statement](#data-availability-statement)
- [Support or Contact](#support-or-contact)

## 1. Molecular dynamics simulation
---

<div style="text-align: justify"> Molecular Dynamics (MD) simulation is a set of molecular simulation methods and a powerful tool for studying condensed state systems. Through molecular dynamics simulation, researchers can get the movement track of the atoms in the system, observe various microscopic details of the atomic movement process, and understand the relationship between the movement of biological macromolecules and their functions, the interaction mechanism between proteins and small molecules. </div>

### 1.1 GROMACS

#### An auto protonation-pdb2gmx run script

<div style="text-align: justify"> The pdb2gmx command of GROMACS when using GROMACS for molecular dynamics simulation, the first command to be used is generally pdb2gmx. This command converts PDB molecular files into GROMACS' unique gro molecular structure file type, and generates molecular topology files at the same time. However, in the process of generating topology files, a protonation process is required, which can define the protonation options of each residue. This step is too cumbersome, so we made this script to auto run protonation for pdb2gmx to prevent dazzle. </div>
<div style="text-align: justify"> <br> </div>
<div style="text-align: justify"> The expect module of Perl is mainly used to practice this function, an example is as follows: </div>

```
#!/usr/bin/perl
use expect;
$exp->expect($timeout,-re=>"Which LYSINE" => sub { $exp->send("1\n"); exp_continue; });      #LYS option
$exp->expect($timeout,-re=>"Which ASPARTIC" => sub { $exp->send("0\n"); exp_continue; });    #ASP option
$exp->expect($timeout,-re=>"Which GLUTAMIC" => sub { $exp->send("0\n"); exp_continue; });    #GLU option
$exp->expect($timeout,-re=>"Which HISTIDINE" => sub { $exp->send("1\n"); exp_continue; });   #HIS option
$exp->interact()
```

The complete automation script is downloaded [here](https://drive.google.com/file/d/1ln_jsnAFGv4qk3abEPBiM5vrZJHcEnf5/view?usp=sharing).
Please contact the author for use permission and the decompression password.

<div style="text-align: justify"> Please note that before executing this script, you need to know the protonation options of each residue. It is recommended to use <a href="https://github.com/mms-fcul/PypKa">pypka method</a> to obtain protonation information. </div>


#### Auto gromacs-result analysis program
<div style="text-align: justify"> An automatic tool for GROMACS based Molecular Dynamics Simulation Analysis.This tool includes two shell scripts. The first script can automatically process the simulation trajectory, including water removal and periodic boundary removal. The second script can automatically export the xvg files of rmsd, rmsf and radius of rotation, as well as the pdb files of short simulation animation. </div>
<div style="text-align: justify"> <br> </div>
<div style="text-align: justify"> The complete automation script is downloaded <a href="https://drive.google.com/file/d/1r_cButINxK7OOXac5bAF70plWuBRDmiq/view?usp=sharing">here</a>. Please contact the author for use permission and the decompression password .</div>

Usage:
```
sh 0auto_traj_process.sh 
sh 1auto_rmsf_analysis.sh
```
<div style="text-align: justify"> Note: You need to manually select the target during runtime o f this tool. If an error occurs, the script will stop automatically .</div>

<div style="text-align: justify"> <br> </div>
<div style="text-align: justify"> There are another  automatic tool <a href="https://heromdanalysis.wordpress.com">HeroMDAnalysis</a> for GROMACS based Molecular Dynamics Simulation Analysis. HeroMDAnalysis is an automagical tool designed to analyze GROMACS based trajectories from molecular dynamics simulations in .xtc format. It can read required .edr , .xtc and .tpr format particle and energy-based coordinate files for biomolecules then perform analysis of various parameters to finally generate suitable high quality images for visualization and publication. </div>

### 1.2 AMBER

#### Online amber tool
<div style="text-align: justify">An online visual amber tool is provided here. Go straight through this <a href="https://cloud.yinfotek.com/">link</a>. The platform uses Amber 20 as the calculation engine, realizing conventional molecular dynamics simulation schemes and rich analysis tools. </div>

#### Script for calculating aMD parameters
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


## 2. Perl
---
<div style="text-align: justify"> Perl is the same as scripting language and doesn't need a compiler and linker to run code. All you have to do is write a program and tell Perl to run it. This means that Perl is ideal for quick solutions to small programming problems and for creating prototypes for large events to test potential solutions. </div>
<div style="text-align: justify"> <br> </div>
<div style="text-align: justify"> Perl has the powerful and flexible characteristics of dynamic language, and also draws on syntax rules from c/c++, basic, Pascal and other languages, thus providing many redundant grammars. In the field of molecular simulation, Perl is often used as a processing script because of its flexibility to text information. The following two examples are provided as a reference for you to learn. </div>

### 2.1 Extract selected residues form protein pdb

<div style="text-align: justify"> In the study of molecular simulation, we may sometimes use the specified partial residues of a protein, so it is very necessary to extract atom of the specified residues. This script mainly solves this problem. This perl script main use the array. Perl array is a list variable that stores scalar values. Variables can be of different types.</div>
<div style="text-align: justify"> <br> </div>
<div style="text-align: justify"> Before using this script, we need to prepare two files, the selected residue file and the pdb file of the protein. The selected residue file is txt format.The sample of selected residue txt file and protein pdb is downloaded <a href="https://drive.google.com/file/d/1Ir7wCGSn9ADX3G9_7Rre0K8wR_CKKQ4g/view?usp=sharing">here</a>.</div>

The following is the usage of this script：
```
perl extr_atom.pl xx.pdb xx.txt
```
Run this script and we will get a XX_ new.pdb, which is all the atom for extracting some residues. The complete script is downloaded [here](https://drive.google.com/file/d/1gfflT5WwTtPfLsbq9Ik9gO0obbwvVqaP/view?usp=sharing).

### 2.2 Extract ligand form complex protein-ligand pdb

<div style="text-align: justify"> This script is fast and convenient tool of extract ligand from the protein complex pdb. This script is usually used in conjunction with <a href="#jump5">Pymol</a> protein prepare script. You can download it from <a href="https://drive.google.com/file/d/1OdRyEdUG_ekzSNBIobFmlqlDb1b8Wsoe/view?usp=sharing">here</a>.</div>

The usage are as follow:
```
perl extr_ligand.pl xx.pdb 
```

NOTE: After running the script, you will get xx_ligand.pdb, you need to manually remove non ligand heteroatoms. 

## 3. Python
---
<div style="text-align: justify"> Python may be one of the few programming languages that can balance simplicity and power at the same time. This is of great benefit to both novices and experts. More importantly, programming in Python is fun. In my research, i often use python to resolve some problem of batch office (including ppt, word, xls and csv),image processing, data visualization and so on. </div>
<div style="text-align: justify"> <br> </div>
<div style="text-align: justify"> Here I will record some of my common scripts and learning experiences for your reference and learning. </div>

NOTE: All scripts of this website are python3 versions.

### 3.1 Python basics 
<div style="text-align: justify"> There is an online website running jupyter notebook, click  <a href="https://colab.research.google.com/notebooks/">here</a>. Input the my library linked (https://github.com/sean28/home) in github block, and you can practice the basics of Python online from here. You can also learn the basic knowledge of Python directly <a href="https://github.com/sean28/home/blob/main/python-basic.ipynb">here</a>.</div>

#### Python reads user input files 

<div style="text-align: justify"> Using python read the user specific files has two way, the first is very simple, for example: </div>


```
#!/usr/bin/python
# -*- coding:utf-8 -*-
import sys
# Get command line parameters
filename = sys.argv[1]   
# Open file
myfile = open(filename,'r') 
# Read each line
for line in myfile.readlines():
# Print each line
        print (line),
# Close file
myfile.close 

```
The usage also simple, commands are as follow:

```
python3 xx.py xx.file
```

<div style="text-align: justify"> The other method is a little complicated, but you can set more function parameters. This method uses a python module argparse. Please refer to <a href="https://docs.python.org/3/library/argparse.html">here</a> for details. The examples are as follows: </div>

```
#!/usr/bin/python
# -*- coding:utf-8 -*-
import argparse
# Create argumentparser object
parser = argparse.ArgumentParser() 
# Add parameters
parser.add_argument('-f', '--filename', type=str, help='input csv filename')
# Analytic parameters
args = parser.parse_args()
# Get parameters
filename = args.filename
print(filename)
```
The usage are as follow:
```
python3 xx.py -f xx.file
```
<div style="text-align: justify"> After mastering the writing method of setting as file input, you can specify Python script to read the specified file later, which is convenient for reading different files, especially suitable for batch processing of a large number of files. </div>

### 3.2 Python processing office 

#### CSV

<div style="text-align: justify"> Comma separated values (CSV, sometimes called character separated values, because the separating character can also be not a comma), its file stores table data (numbers and text) in plain text. Plain text means that the file is a sequence of characters without data that must be interpreted like binary numbers. CSV file is composed of any number of records, which are separated by some line break; Each record consists of fields. The separator between fields is other characters or strings, and the most common is comma or tab. Generally, all records have exactly the same sequence of fields. They are usually plain text files.</div>
<div style="text-align: justify"> <br> </div>
<div style="text-align: justify">For the field of bioinformatics, it is very convenient to use the plain text format of csv instead of Excel to process some data, and there is no need to import additional modules or expansion packages that support reading xls and xlsx.</div>

(1) Fast delete data of blank script

We can use this script fast delete the data with blank column in the csv table. You can download it from <a href="https://drive.google.com/file/d/1bfIP5UpzTOtHbYDpVnugB09Mb3q2Dm4Z/view?usp=sharing">here</a>.

The usage are as follow:
```
python3 delete_csv_black.py -f xx.csv
```
(2) Batch fetch pdb and ligand
<br>
to be continue...

#### EXCEl

to be continue...

#### PPT 
(1) Auto delete ppt
<br>
This is a sample script that can quickly delete the first and last pages of all ppt files in the current directory. You can download it from <a href="https://drive.google.com/file/d/1Ka4PSZs7jii4buE0355zzSHq46uL4hdy/view?usp=sharing">here</a>.

The usage are as follow:
```
python3 delete-ppt.py
```

(2) Merge ppt
<br>
This is a script that can quickly merge all ppts under the current folder. You can download it from <a href="https://drive.google.com/file/d/1Yt9MWSQzN0e1J2xKztOANS2pZ4gupBsr/view?usp=sharing">here</a>.
```
python3 merge-ppt.py
```



### 3.3 Python processing image 
<div style="text-align: justify"> This part of the script is cumbersome and does not include the function of user specified file input. It is temporarily presented in code for reference only. </div>

#### Vertical merge png

```
#!/usr/bin/python
# -*- coding:utf-8 -*-
# Created By Sean -- MUST ;
from PIL import Image

img1 = Image.open( "./image1.png")
img2 = Image.open( "./image2.png")

img3 = Image.new('RGB', (img1.size[0], img1.size[1] + img2.size[1])
img3.paste(img1, (0, 0))
img3.paste(img2, (0, img1.size[1]))
img3.paste(img3, (0, img1.size[1] + img2.size[1]))

img3.save("./result_merge_ver.png")
```

#### Horizontal merge png

```
#!/usr/bin/python
# -*- coding:utf-8 -*-
# Created By Sean -- MUST ;
from PIL import Image

img1 = Image.open( "./image1.png")
img2 = Image.open( "./image2.png")

img3 = Image.new('RGB', (img1.width + img2.width, img1.height))
img3.paste(img1, (0, 0))
img3.paste(img2, (img1.width, 0))
img3.save("./result_merge_hor.png")
```


## 4. R
---

### 4.1 Parallel kmeans scripts 
<div style="text-align: justify">R is a powerful scripting language for mapping and data visualization, which can execute a large number of mathematical models and algorithms. However, due to its low system execution efficiency, it will be difficult to deal with the problem of large amount of data. Here is a case of parallel kmeans clustering for everyone to learn. </div>
<div style="text-align: justify"> <br> </div>
<div style="text-align: justify">For the parallel running script of kmean mean value, the best of the 10 running results is better selected. After testing, it does not affect the operation results, greatly speeds up the operation speed and saves the script operation time. It is recommended to use when calculating large data sets.  </div>
<div style="text-align: justify"> <br> </div>
<details>
<summary> Click to view code </summary>
<pre><code>
#####By Sean from MUST#####
#Induce parallel package
library(parallel)

#Define number of cpu core (default: total -2)
nw <- detectCores()-2
cl <- makeCluster(nw)

#Define nstart 
nstart <- 10
nstartv <- rep(ceiling(nstart / nw), nw)

#read lig_noh_pos.txt
data <- read.table("data.txt")

#run clusterApply
data_km <- clusterApply(cl, nstartv,
        function(n, x) kmeans(x, 1000, nstart=n, iter.max=100),
        data)
        
#Pick the best result
i <- sapply(data_km , function(data_km) data_km $tot.withinss)
data_km  <- data_km [[which.min(i)]]
print(data_km$tot.withinss)
per_atom_rmsd<-sqrt((data_km$withinss/(data_km$size-1))/2914) 
summary(data_km$size)
summary(per_atom_rmsd)
</code></pre>
</details>

### 4.2 ggplot2 
<a href="url"><img src="https://ggplot2.tidyverse.org/logo.png" align="center" height="56" ></a>
<div style="text-align: justify"> Ggplot is an R software package used to draw statistical graphs. It is an important tool to visualize data analysis, supported by a set of syntax behind it. The core idea of ggplot2 is to separate drawing and data, and separate data related drawing from data independent drawing. </div>
<div style="text-align: justify"> <br> </div>
<div style="text-align: justify"> In the first, let us learn about the installation of ggplot2, as shown below: </div>

Installation：
```
# Check the installed package of R
library()
# install ggplot2
install.packages("ggplot2")
# import ggplot2
library(ggplot2)
```
<div style="text-align: justify"> Basic usage of ggplot2 </div>
<div style="text-align: justify"> <br> </div>
to be continue...

### 4.3 Binder 
<div style="text-align: justify"> There is a online tool can run the R, binder. The <a href="https://mybinder.org/">binder</a> can directly configure the environment of GitHub as a docker image, and then start it in the cloud. With Binder,we can open those notebooks in an executable environment. I have deployed R in my public repository, and through this <a href="https://mybinder.org/v2/gh/sean28/home/HEAD">link</a>, you can learn and practice the basic knowledge of R language or python online.</div>
<div style="text-align: justify"> <br> </div>

### 4.4 TRAPP Multiple Comparison Script 
<div style="text-align: justify"> TRAnsient Pockets in Proteins (TRAPP) is a tool that allows the exploration of different protein conformations, the analysis of binding pocket flexibility and dynamics, and the extraction of spatial and physicochemical information on the binding pocket confor-mations (J Chem Inf Model. 2020 Mar 23;60(3):1685-1699). Through this <a href="https://trapp.h-its.org/trapp">link</a>, you can learn and use the TRAPP webserver online.</div>
<div style="text-align: justify"> <br> </div>
<div style="text-align: justify"> The analysis chart derived from TRAPP only shows the scoring of one system, and cannot compare multiple groups of systems. Now we have developed an extended tool to realize the comparative analysis of multiple systems. This program need to extract data from the TRAPP analysis results, and then use this script for analysis. Now we only provide two groups of system comparative analysis tools. </div>
![This is an image](/sample.png)
<div style="text-align: justify"> This tool is now open source and you can refer to it from <a href="https://github.com/sean28/TRAPP-Multiple-comparison.git">here</a>. Download this tool from this <a href="https://github.com/sean28/TRAPP-Multiple-comparison/archive/refs/heads/main.zip">link</a>. If you need to use it, please indicate the source. Articles using this tool have been published, please refer to this link (https://pubmed.ncbi.nlm.nih.gov/36232570/).</div>
<div style="text-align: justify"> <br> </div>

## 5. Pymol 

<div style="text-align: justify"> PyMOL is a molecular 3D structure display software, which is suitable for creating high-quality 3D structure images of small molecules or biological macromolecules. The content of this issue is to share the quick reference manual of common commands in the use of PyMOL. If you need the original version, you can browse it on the official website <a href="https://pymolwiki.org/index.php/Main_Page">(pymolwiki)</a>.</div>


|Description           |        Command              |
|          :---        |           :---              |
|Start record operation|  log_open xx.pml|
|Close record operation|  log_close      | 
|remove solvent        | remove solvent  |
|Select Na+            | select selectionname, r. na\\+|
|Select Cl-            | select selectionname, r. cl\\-|
|Delete selection      |  delete ,selection|
|Merge selection       |  select ,new_selection, selection1 + selection2  | 
|Structural superposition | cealign, align, super, pairfit, fit  selection, target  <a href="https://github.com/sean28/home/blob/ba4be1a571124c8a258852e75c837058bdf112d8/Superpositions.pdf">(Detailed description)</a>|
|Label emerges         | set label_position, (0, 0, 40) |
|Print label           |  iterate  selection and name CA, print (label) |
|Print resi and resn   |  iterate  selection and name CA, print (resi, resn)  | 
|Rotating angle of view| turn x, 180 or turn y, 180|
|Select residues around 5 angstroms of ligand | select 5A_resi, byres selection_ligand_name around 5|
|Estimate the x, y and z of protein  |  get_extent <selection_name>|
|Move the screen 10A along the X axis ̊   |  move x, 10  | 
|View the number of residues  |  print len(cmd.get_model("poly").get_residues()) | 
|View the number of atoms     |  count_atoms <selection_name>| 
|Altering atom coordinates    |  alter_state 1, <selection_name>, y=y+10.0 |
|Cacluate RMSD                |  rms (selection), (target-selection)|
|select anything shown as a lines/dots/spheres |  select rep lines/dots/spheres|

 
There are also some simple scripts that can easily implement some functions, and  You can download them here for free.

(1) [Obtain protein residue sequence information. ](https://github.com/sean28/home/blob/bddec2088f8db6c096d0f62b22bc9d30f810eeff/get_fasta.pml)

(2) [Hydrogenation of protein.](https://github.com/sean28/home/blob/bddec2088f8db6c096d0f62b22bc9d30f810eeff/add_H.pml)

(3) [Remove solvent molecules.](https://github.com/sean28/home/blob/bddec2088f8db6c096d0f62b22bc9d30f810eeff/remove_water.pml)

to be continue...



## 6. Free Energy 
---
<div style="text-align: justify"> In chemistry, the lower the free energy is, the greater the affinity between the receptor and the ligand is, and the more likely the molecular docking reaction is to occur. That is, the lower the binding free energy, the easier the key and lock will be stuck together, and the more effective the drug will be. It can destroy the normal function of protein more effectively. The prediction of binding free energy that is meaningful for the optimization of lead compounds needs to be within 1kcal/mol (~0.04% of the total energy).</div>
<div style="text-align: justify"> <br> </div>
<div style="text-align: justify"> Generally speaking, the binding affinity of drugs is between -8.28kcal/mol and -12.43kcal/mol. For every 1 kcal/mol difference in binding free energy, the activity difference is about 6 times. A 10 fold change in affinity is equivalent to a change in the binding free energy of 1.4kcal/mol.</div>
<div style="text-align: justify"> <br> </div>
Molar concentration unit table:

<p align="left">
 <img src="/Molar-unit-table.png" width="500" >
</p>


<div style="text-align: justify">In addition, there are calculation tools of free energy is  <a href="https://drive.google.com/file/d/1x8zNoy30bsR6UmWtqnQAziPVs-8cLuHL/view?usp=sharing">here</a>. For conversion between Binding free energy ΔG and affinity (IC50/kd/ki). </div>

The formula is: ΔGbinding = RT·ln Kdissociated = RT·lnKd ≈ RT·lnIC50 = −RT·pIC50

Usage:

Ic50/ki/kd to ΔGbinding:
```
python cal_gbinding.py
Please enter the temperature (K)：
Please enter the IC50 (μM)：
```

ΔGbinding to Ic50/ki/kd:
```
python cal_ic50.py
Please enter the temperature (K)：
Please enter the Gbinding (kcal/mol)：
```
## 7. Others
---
### 7.1 Automatically building scientific research environment

<div style="text-align: justify"> For many novices who are just beginning in bioinformatics, it is a headache to explore the installation environment of different software. CONDA can create different virtual environments and install them into different virtual environments according to different software, which is not easy to conflicts due to the dependencies of different programs. It only needs to call different virtual CONDA environments to call different software. This method can use in linux and mac os.</div>
<div style="text-align: justify"> <br> </div>
Create a new virtual environment under the specified directory, and enter the command:

```
# Python 
conda create --prefix=./conda_work/envs/python3 python=3.8
conda create --prefix=./conda_work/envs/python2 python=2.7

# amber 
conda create --prefix=./conda_work/envs/amber  --no-default-packages
conda activate amber
conda install -c conda-forge ambertools

# gromacs
conda create --prefix=./conda_work/envs/gromacs --no-default-package
conda activate gromacs
conda install -c bioconda gromacs

# R
conda create -n R3.5
source activate R3.5
conda install r-base=3.5.1
## R packages usually need to start with r-
conda install r-ggplot2
## If it cannot be found, you can use this command to search the corresponding R package. Anaconda here is the path of the original CONDA, not in r3.5 environment
anaconda search -t conda r-ggplot2
## Show Chanel of this package
anaconda show BioBuilds/r-ggplot2
## Install according to Anaconda show
conda install --channel https://conda.anaconda.org/BioBuilds r-ggplot2

# Openbabel (Note: Environment requiring Python 2.7)
conda install -c openbabel openbabel

```

Use the command to view the currently owned virtual environment:
```
conda info --envs
```

Set the shortcut opening method for each virtual environment and write the following lines into the environment file (Linux: ~/.bashrc MACOS: ~/.zshrc).

```
alias work="conda activate ../python3"
alias py2="conda activate ../python2"
alias amber="conda activate ../amber"
alias gromacs="conda activate ../gromacs"
```
Note: The path (../python3) should be a full path.

Delete virtual environment:

Input for example:
```
conda remove --name en_name --all
```
<div style="text-align: justify"> In the last, it is worth noting that if you want to use multiple software in the same virtual environment, and the software will call each other, you need to deploy multiple software in the same environment. For example, when using Jupiter notebook, Python and R, other situations will not be listed one by one.</div>
<div style="text-align: justify"> <br> </div>

### 7.2 Online calculator
<div style="text-align: justify"> When doing computational chemistry, the conversion of various units is a headache, especially when accurate values are really needed. The online calculator developed by jerkwin is very easy to use. His github link is here (https://github.com/Jerkwin/gmxtools). However, the webpage of the calculator crashed and could not be logged in. So, I redeployed this online tool in my repository and you can use this tool for free by clicking <a href="https://sean28.github.io/Online_calculator/">here</a>. Besides, I added a new unit conversion in this calculator, nM(nmol/L). If there is infringement, please contact to delete.</div>
<div style="text-align: justify"> <br> </div>

### 7.3 File format converter
<div style="text-align: justify">The conversion between various chemical formats is a headache question. Different chemical formats have different uses.For example, In the field of molecular dynamics simulation, the format is used as follows:</div>


|Name|Usage|
|:---|:---|
|Topology file.     |The coordinate file records the three-dimensional coordinates of all atoms in the simulation system. Format include: amber (prmtop, pram7) gromacs (top) |
|Coordinate file.   |Topological files record the connection relations of atoms and molecular mechanical parameters. Format include: amber (inpcrd) gromacs (gro)|
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

### 7.4 Online chemical structure drawing and editing tool
<div style="text-align: justify">ChemDraw JS is a chemical structure drawing and editing tool designed to help you create high-quality chemical drawings. ChemDraw JS consists of a drawing toolbar and document window. The document window lets you draw and edit chemical structures using the various drawing tools available in the ChemDraw JS toolbar. This section introduces the ChemDraw JS user interface and ChemDraw JS toolbar.</div>
<div style="text-align: justify"> <br> </div>
<div style="text-align: justify">Click <a href="https://chemdrawdirect.perkinelmer.cloud/js/sample/index.html#">here</a> to jump to this tool.</div>

## Data Availability Statement

If you use the data of this website, please indicate the source of the reprinted article. 

## Support or Contact

<div style="text-align: center"> Team of drug design and discovery </div>
<div style="text-align: center"> Dr. Neher’s Biophysics Laboratory for Innovative Drug Discovery </div>
<div style="text-align: center"> Macau University of Science and Technology </div>
<div style="text-align: center"> Avenida WaiLong, Taipa, Macau(SAR) </div>
<div style="text-align: center"> <a href="mailto:xjyao@must.edu.mo">Professor Dr. Yao Xiao-jun </a> </div>
<div style="text-align: center"> <a href="mailto:sean28299@gmail.com">Dr. Sean</a> </div>
<div style="text-align: justify"> <br> </div>

<div align=center><img src="/Batman.png"></div>
---
