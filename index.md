# 🧬 Sean Drug Design and Discovery

> **Welcome to the official site of Dr. Sean's group.**  
> We focus on *computer-aided drug design (CADD)*, *molecular simulation*, *computational chemistry*, and *bioinformatics*. This site serves as a collection of learning notes and self-developed tools.

<p align="center">
  <img src="/cadd.png" alt="CADD banner" width="60%">
</p>

---

## 📌 Introduction

This repository is maintained by **Dr. Sean**, aiming to provide practical scripts, toolkits, and guides for molecular modeling and data analysis. Most content is developed in real research projects and validated by peer-reviewed publications.

If any content infringes, please contact for deletion.  
Contact: [sean28299@gmail.com](mailto:sean28299@gmail.com)

---

## 🗂️ Table of Contents

### 🧪 1. Molecular Dynamics Simulation
- **[1.1 GROMACS](#11-gromacs)**  
  - An auto protonation-pdb2gmx run script  
  - Auto gromacs-result analysis program  
- **[1.2 AMBER](#12-amber)**  
  - Online amber tool  
  - Script for calculating aMD parameters

### 🧬 2. PDB Operations
- Extract selected residues  
- Extract ligand from protein-ligand complex  
- Fetch PDB ID validation report

### 🔬 3. Pymol
- Common commands and handy scripts

### 📉 4. Free Energy Calculations
- ΔG and IC50 interconversion  
- Molar unit reference table

### 🐍 5. Python Programming
- File I/O, argparse, csv/xls/ppt/image  
- Scripts for bioinformatics and automation

### 📊 6. R Programming
- Parallel kmeans  
- ggplot2 & TRAPP multiple system comparison

### 🧰 7. Other Tools
- Auto environment build with Conda  
- Online calculator & color map tool  
- File converter, ChemDraw JS, etc.

---

## 📤 Data & Usage Statement

Please cite or acknowledge when using data/scripts from this site. All downloadable scripts are provided **as-is** for academic use.

---

## 📬 Contact

<p align="center">
  <strong>Team of Drug Design and Discovery</strong>  
  📧 Prof. Yao Xiao-jun — <a href="mailto:xjyao@must.edu.mo">xjyao@must.edu.mo</a>  
  📧 Dr. Sean — <a href="mailto:sean28299@gmail.com">sean28299@gmail.com</a>  
</p>

<p align="center">
  <img src="/Batman.png" alt="Sean logo" width="30%">
</p>

---

## 🧪 1. Molecular Dynamics Simulation

Molecular Dynamics (MD) simulation is a powerful method to explore atomic-level movements and understand the functional mechanisms of biomolecular systems.

---

### ⚙️ 1.1 GROMACS

#### 🧼 Auto Protonation & `pdb2gmx` Run Script

`pdb2gmx` is often the first command used in GROMACS, responsible for converting a `.pdb` file to `.gro` format and generating topology files. The residue protonation options required can be tedious, so we provide an automated script using Perl's `Expect` module.

```perl
#!/usr/bin/perl
use expect;
$exp->expect($timeout,-re=>"Which LYSINE" => sub { $exp->send("1\n"); exp_continue; });
$exp->expect($timeout,-re=>"Which ASPARTIC" => sub { $exp->send("0\n"); exp_continue; });
$exp->expect($timeout,-re=>"Which GLUTAMIC" => sub { $exp->send("0\n"); exp_continue; });
$exp->expect($timeout,-re=>"Which HISTIDINE" => sub { $exp->send("1\n"); exp_continue; });
$exp->interact()
```

> 📥 **[Download Script](https://drive.google.com/file/d/1ln_jsnAFGv4qk3abEPBiM5vrZJHcEnf5/view?usp=sharing)**  
> 🔐 Contact the author for password & permission.  
> 🔍 Protonation recommendation: [PypKa tool](https://github.com/mms-fcul/PypKa)

---

#### 📊 Auto GROMACS Result Analysis

We provide two shell scripts for post-MD trajectory processing:

- 🧹 `0auto_traj_process.sh` — water removal, periodic boundary fix  
- 📈 `1auto_rmsf_analysis.sh` — generate RMSD/RMSF/Rg `.xvg` and short trajectory `.pdb`

```bash
sh 0auto_traj_process.sh
sh 1auto_rmsf_analysis.sh
```

> 📥 **[Download Scripts](https://drive.google.com/file/d/1r_cButINxK7OOXac5bAF70plWuBRDmiq/view?usp=sharing)**  
> ❗ Requires manual selection during run. Script will auto-stop if error occurs.

> 🧠 For a more visual solution, check: [HeroMDAnalysis](https://heromdanalysis.wordpress.com)

---

### 🧬 1.2 AMBER

#### 🌐 Online Amber Tool

A web-based AMBER20 platform is available via:  
🔗 [https://cloud.yinfotek.com](https://cloud.yinfotek.com)

This online platform supports conventional MD simulation workflows and rich result analysis.

---

#### ⚡ aMD Parameter Calculation Script

Accelerated MD (aMD) modifies potential energy to enhance sampling. You’ll need to calculate:

| Parameter   | Description                        | Formula |
|-------------|------------------------------------|---------|
| EthreshP    | Total potential energy threshold   | `E(tot) = EPtot + 0.16 × atom_num` |
| alphaP      | Boost factor for total potential   | `alphaP = 0.16 × atom_num` |
| EthreshD    | Dihedral energy threshold          | `E(dih) = DIHED + 4 × resi_num` |
| alphaD      | Boost factor for dihedral energy   | `alphaD = 0.2 × 4 × resi_num` |

Python script:

```python
#!/usr/bin/python
# Calculate EthreshP / alphaP / EthreshD / alphaD
EPtot = float(input("EPtot (kcal/mol): "))
DIHED = float(input("DIHED (kcal/mol): "))
atom_num = int(input("Number of atoms: "))
resi_num = int(input("Number of solute residues: "))

EthreshP = round(EPtot + 0.16 * atom_num, 2)
alphaP = round(0.16 * atom_num, 2)
EthreshD = round(DIHED + 4 * resi_num, 2)
alphaD = round(0.2 * 4 * resi_num, 2)

print("EthreshP =", EthreshP)
print("alphaP =", alphaP)
print("EthreshD =", EthreshD)
print("alphaD =", alphaD)
```

---

## 🧬 2. PDB Operations

PDB (Protein Data Bank) files store 3D structural information of proteins, nucleic acids, and other biomolecules. This section contains handy scripts for efficiently manipulating `.pdb` files.

---

### 🔎 2.1 Extract Selected Residues

Sometimes only specific residues from a protein are needed for simulations. This Perl script extracts atoms of target residues.

```bash
perl extr_atom.pl input.pdb residue_list.txt
```

- 📂 Sample input files: [Download](https://drive.google.com/file/d/1Ir7wCGSn9ADX3G9_7Rre0K8wR_CKKQ4g/view?usp=sharing)
- 📥 Script: [Download extr_atom.pl](https://drive.google.com/file/d/1gfflT5WwTtPfLsbq9Ik9gO0obbwvVqaP/view?usp=sharing)
- 📌 Output: `input_new.pdb`

---

### 💊 2.2 Extract Ligand from Complex PDB

Quickly extract ligand atoms from protein-ligand complexes. Useful for ligand-only preparation.

```bash
perl extr_ligand.pl complex.pdb
```

> 📥 [Download Script](https://drive.google.com/file/d/1OdRyEdUG_ekzSNBIobFmlqlDb1b8Wsoe/view?usp=sharing)  
> ⚠️ Output: `complex_ligand.pdb` (remove unwanted HETATM manually)

---

### 📑 2.3 Fetch PDB ID Validation Report

Python script for batch downloading validation PDFs by entering multiple PDB IDs.

- Prompt: `Enter comma-separated PDB IDs (e.g., 6LU7,1CBS,2PTC)`
- Output folder: `validation_reports/`

> 📥 [Download Script](https://drive.google.com/file/d/1fyBODIrKMvWLFQbuy1XsYgsvHXO1VC7b/view?usp=sharing)`

