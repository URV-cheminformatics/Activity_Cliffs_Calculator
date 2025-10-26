# 🧮 Activity Cliff Calculator

## A Python tool to identify and analyze **activity cliffs** — pairs of structurally similar compounds with large potency differences.

It takes as input an .xls or .xlsx file containing molecular data and automatically detects activity cliffs based on the differences in activity (pIC50) and chemical similarity (Tanimoto coefficient). It computes **Tanimoto similarity**, **ΔpIC50**, and **disparity** using multiple RDKit fingerprints, and supports **filtering, Excel export with molecule images**, and **inline SVG visualization** directly in Jupyter or VS Code.

---

## 📘 Overview

This tool detects **activity cliffs**, defined as compound pairs (or groups) that share high structural similarity but exhibit large differences in potency.  
It automates the generation of molecular pairs, similarity computation, filtering, visualization, and export.

---

📂 Input requirements

The input file must contain at least the following columns:
| Column | Description |
|:--------|:-------------|
| `standardize_smiles` | Canonical SMILES string (used to compute fingerprints) |
| `PubMedID` | Identifier used to group compounds from the same publication |
| `pIC50` | Biological activity value (−log₁₀(IC₅₀)) |

Additional columns (e.g. `Compound`, `Target`, etc.) are preserved in the output.

---

## ⚙️ Main features

- 🧬 Supports **7 RDKit fingerprint types**:`morgan`, `feature_morgan`, `rdk`, `maccs`, `pattern`, `atompair`, `torsion`
- Fully compatible with **RDKit 2025.01+**, but adapts the fingerprint generation for several RDKit versions
- 🔁 Pair generation grouped by `PubMedID` (or other columns)
- 📊 Flexible filtering by **Tanimoto**, **ΔpIC50**, and **Disparity**
- 💾 Excel export with **embedded molecule images**
- 🧠 Inline molecule preview (SVG) in VS Code / Jupyter
### 🔁 Pair Generation
- Generates all compound pairs **within the same group** (default: `PubMedID`)
- Computes:
  - **Tanimoto similarity**
  - **Activity difference (ΔpIC50)**
  - **Disparity** = {Delta pIC50}/(1 - Tanimoto)
  - Handles identical molecules (`Tanimoto = 1`) gracefully

 ### 🎚️ Filtering Options
You can filter the resulting pairs using any combination of:
- `tanimoto_min` → minimum similarity
- `activity_diff_min` → minimum ΔpIC50
- `disparity_min` → minimum disparity value  

Each filter is **optional** — if set to `None`, it is ignored.

---

## 📊 Visualization and Export

### 📈 Interactive Plot
- Creates an **interactive scatter plot** (`ΔpIC50` vs. `Tanimoto`) using Plotly  
- One plot per selected fingerprint type

### 💾 Excel Export
- Exports all results to `.xlsx` format with **molecule images embedded** (PNG)
- Automatically adjusts **row height** and **column width** to match image size  
- Configurable image size via `image_size` parameter (e.g., 75, 100, 150 px)

### 🧠 Inline Visualization
- `show_molecule_table()` displays molecules as **inline SVGs** directly in VS Code or Jupyter Notebook
- Automatically limits drawing to the first *n* rows (default: 20) for performance  
- Adjustable preview size with `img_size=(width, height)`

---

## 🧰 Example Workflow

```python
from activity_cliffs_utils import (
    fp_as_bitvect, generate_pairs,
    export_activity_cliffs_to_excel, show_molecule_table, mol_to_image_bytes, smiles_to_svg
)
from rdkit import Chem
import pandas as pd
import plotly.express as px


# 1️⃣ Load input data
df = pd.read_excel("M-pro_Inhibitors.xls")

# 2️⃣ Compute fingerprints (feature_morgan by default)
df["fp_feature_morgan"] = df["standardize_smiles"].apply(
    lambda s: fp_as_bitvect(Chem.MolFromSmiles(s)) if Chem.MolFromSmiles(s) else None
)

# 3️⃣ Generate compound pairs (grouped by PubMedID)
result_df = generate_pairs(df, group_col="PubMedID")

# 4️⃣  Create the interactive chart
fig = px.scatter(
    result_df,
    x='Tanimoto',
    y='pIC50_diff',
    title='Activity difference vs Tanimoto Similarity',
    labels={'Tanimoto': f'Tanimoto Similarity (feature_morgan)', 'pIC50_diff': 'Activity difference (ΔpIC50)'},
    hover_data=['Compound1', 'Compound2', 'Disparity']
)
# Adjust the chart dimensions
fig.update_layout(
    height=700  # You can change this value to adjust the height
)
# Show the graph
fig.show()
```

<img width="1733" height="472" alt="imatge" src="https://github.com/user-attachments/assets/517956f7-8cd7-43a9-9615-3b4df3555046" />

```
# 5️⃣ Apply filters
tanimoto_min = 0.8
activity_diff_min = 1.5
mask = (result_df["Tanimoto"] >= tanimoto_min) & (result_df["pIC50_diff"] >= activity_diff_min)
filtered_df = result_df[mask].copy()

# 6️⃣  Add molecule images for Excel export and preview
filtered_df["Mol1_img"] = filtered_df["SMILES1"].apply(lambda s: mol_to_image_bytes(s, (150,150)))
filtered_df["Mol2_img"] = filtered_df["SMILES2"].apply(lambda s: mol_to_image_bytes(s, (150,150)))
filtered_df["Mol1_svg"] = filtered_df["SMILES1"].apply(lambda s: smiles_to_svg(s, size=(120,120)))
filtered_df["Mol2_svg"] = filtered_df["SMILES2"].apply(lambda s: smiles_to_svg(s, size=(120,120)))

# 7️⃣  Export to Excel
export_activity_cliffs_to_excel(filtered_df, "disparity_results_fp_morgan.xlsx", image_size=150)

# 8️⃣ Display molecule preview
show_molecule_table(filtered_df, max_rows=30, img_size=(120,120))
```
<img width="1153" height="354" alt="imatge" src="https://github.com/user-attachments/assets/240c8655-6adf-4f0f-9d86-890e5b270ce8" />


🧾 Output Columns
| **Column** | **Description** |
|:------------|:----------------|
| `PubMedID` | Group identifier |
| `Compound1`, `Compound2` | Compound names |
| `pIC50_1`, `pIC50_2` | Activity values |
| `pIC50_diff` | Absolute difference in activity |
| `Tanimoto` | Structural similarity |
| `Disparity` | Activity cliff magnitude |
| `Fingerprint` | Fingerprint type used |
| `Mol1_img`, `Mol2_img` | Molecule images (PNG, Excel export only) |

---

🧩 Highlights

✅ Compatible with RDKit 2025+

✅ Supports 7 fingerprint types

✅ Flexible filtering and visualization

✅ Integrated Excel export with molecule images

✅ Inline SVG previews for quick inspection

✅ Clear, English docstrings and maintainable code

---

## 🧑‍💻 Authors and License

Developed by Santi Garcia-Vallvé. Universitat Rovira i Virgili. 
Department of Biochemistry and Biotechnology.
Cheminformatics and Nutrition (QiN) Research Group

Parts of the code structure and documentation were drafted with assistance from ChatGPT (OpenAI), under human supervision and subsequent review.

Licensed under the MIT License.
