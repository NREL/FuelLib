import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os
import re
import sys

# Add the FuelLib directory to the Python path
fuellib_dir = os.path.dirname(os.path.dirname(__file__))
sys.path.append(fuellib_dir)
import GroupContributionMethod as gcm

fuel_name = "posf10325"

fuel = gcm.groupContribution(fuel_name)

# Read gcxgc file
gcxgcFile = os.path.join(fuel.gcxgcDir, f"{fuel.name}_init.csv")
df = pd.read_csv(
    gcxgcFile,
)

# Get column names
colNames = df.columns.tolist()

# Classify rows based on the first column:
# if df[colNames[0]] contains Toluene, Benzene or aromatic, then classify as aromatic
aromatic = df[colNames[0]].str.contains(
    "Toluene|Benzene|Aromatic", case=False, na=False
)
# if df[colNames[0]] contains n-C, then classify as n-alkane
n_alkane = df[colNames[0]].str.contains("n-C", case=False, na=False)
# if df[colNames[0]] contains isoparaffin, then classify as iso-alkane
isoalkane = df[colNames[0]].str.contains("Isoparaffin", case=False, na=False)
# if df[colNames[0]] contains cycloparaffin, then classify as cyclo-alkane
cycloalkane = df[colNames[0]].str.contains("Cycloparaffin", case=False, na=False)

# Append classification as a new column
family_names = ["n-alkane", "iso-alkane", "cyclo-alkane", "aromatic"]
df["Family"] = np.select(
    [n_alkane, isoalkane, cycloalkane, aromatic], family_names, default="unknown"
)


# Determine carbon number by row:
def determine_carbon_number(compound):
    if "Toluene" in compound:
        return 7
    elif "benzene" in compound.lower():
        # Extract the number after "C" and add 7
        match = re.search(r"C(\d+)", compound)
        if match:
            try:
                carbon_number = int(match.group(1))
                return carbon_number + 6
            except ValueError:
                return np.nan  # Handle cases where extraction fails
        else:
            return np.nan  # No match found
    else:
        # Extract the number after "C" in either format
        match = re.search(r"C(\d+)", compound)
        if match:
            try:
                return int(match.group(1))
            except ValueError:
                return np.nan  # Handle cases where extraction fails
        else:
            return np.nan  # No match found


# Apply the function to the column and append as a new column
df["nC"] = df[colNames[0]].apply(determine_carbon_number)

# Remove rows <= 0.01 in weight % column at max(nC)
df = df[df[colNames[1]] > 0.01]

# Plotting parameters
spacing = [-0.2985, -0.099, 0.099, 0.2985]
colors = {
    "n-alkane": "#063C61",
    "iso-alkane": "#2980B9",
    "cyclo-alkane": "#91BCD8",
    "aromatic": "#7f7f7f",
}

# Bar plot of the number of compounds in each family
plt.figure(figsize=(7, 5))
N = df.nC.unique()
for k, family in enumerate(family_names):
    nC = df[df["Family"] == family].nC
    weight = df[df["Family"] == family][colNames[1]]

    # check duplicate nC values
    if len(nC) != len(set(nC)):
        # If there are duplicates, sum the weights for each nC
        df_grouped = (
            df[df["Family"] == family].groupby("nC")[colNames[1]].sum().reset_index()
        )
        nC = df_grouped["nC"]
        weight = df_grouped[colNames[1]]
    plt.bar(
        nC + spacing[k], weight, label=family, alpha=1, color=colors[family], width=0.2
    )
    plt.xticks(df.nC.unique(), fontsize=14)
    plt.xlim(min(N) - 0.5, max(N) + 0.5)

plt.xlabel("Carbon Number", fontsize=16)
plt.xticks(df.nC.unique(), fontsize=14)
plt.yticks(fontsize=14)
plt.ylabel("Weight %", fontsize=16)
plt.title("Fuel Composition", fontsize=16, fontweight="bold")
plt.legend(fontsize=14)
plt.tight_layout()
plt.show()
