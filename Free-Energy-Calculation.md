# Free Energy 
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
