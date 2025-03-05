Calculate activity cliffs between chemical compounds.

Compute activity Cliffs (pairs or groups of structurally similar compounds that are active against the same target but have large differences in potency) from an xls file containing a set of molecules. The columns "standardize_smiles", "PubMedID" and "pIC50" are mandatory. 
The program uses rdkit to compute Morgan fingerprints of each molecule. When the PubMedID is the same, for each pair of molecules a disparity value is calculated as disparity = pIC50_diff / (1 - tanimoto).

The program needs pandas, numpy, tqdm, rdkit, combinations and plotly.express libraries

Example on an input file https://www.mdpi.com/article/10.3390/ijms23010259/s1 from Macip G, Garcia-Segura P, Mestres-Truyol J, Saldivar-Espinoza B, Pujadas G, Garcia-Vallvé S. A Review of the Current Landscape of SARS-CoV-2 Main Protease Inhibitors: Have We Hit the Bullseye Yet? Int J Mol Sci. 2021 Dec 27;23(1):259. doi: 10.3390/ijms23010259# Activity_Cliffs_Calculator

![newplot](https://github.com/user-attachments/assets/613bf395-c978-4fe0-a651-1ff75312ec26)
