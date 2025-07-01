import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os
import GroupContributionMethod as gcm

# -----------------------------------------------------------------------------
# Calculate mixture properties from the group contribution properties
# -----------------------------------------------------------------------------

# HEFA fuels from various feedstocks (see fuelData/propertiesData for fuels)
fuel_names = ["hefa-mfat", "hefa-came", "hefa-tall"]
conv_fuel_name = "jet-a"
blends = np.linspace(0,100,100) # Weight percentages of hefa in blend

# Properties to plot
prop_names = ["Density","Viscosity"]

# Plotting parameters
fsize = 18
ticksize = 18
line_thickness = 4
marker_size = 75

# droplet specs
drop = {}
drop["d_0"] = 100 * 1e-6  # initial droplet diameter (m), note: size doesn't matter
drop["r_0"] = drop["d_0"] / 2.0  # initial droplet radius (m)

# Line specifications for plotting
def linespecs(name):
    if "came" in name:
        return "#7f7f7f", "o"  # 50% Gray
    elif "mfat" in name:
        return "#333333", "o"  # Dark Gray
    elif "tall" in name:
        return "#2980B9", "o"  # Primary Blue
    else:
        return "#2980B9", "o"  # Fallback: Primary Blue


# Legend labels for plotting
def leglab(name):
    if "hefa" in name:
        return name.upper()
    else:
        return "Jet-A"

# y-axis label
def ylab(prop_name, temp):
    if prop_name == "Density":
        return rf"{prop_name} at {temp:g} °C [g/cm$^3$]"
    elif prop_name == "Viscosity":
        return rf"{prop_name} at {temp:g} °C [mm$^2$/s]"

def getPredAndData(fuel_name, prop_name, blend):
    blend = np.array(blend)*1e-2  # Convert to weight percent

    # Get the fuel properties based on the GCM
    fuel = gcm.groupContribution(fuel_name, "hefa")
    jetA = gcm.groupContribution(conv_fuel_name)

    data_file = "hefa-jet-a-blends.csv"
    dataPath = os.path.join(fuel.fuelDataDir, "propertiesData")
    data = pd.read_csv(os.path.join(dataPath, data_file), skiprows=[1])
    col = f"{prop_name}_{fuel_name[5:].upper()}"
    prop_data = data[col]
    blend_data = data["HEFA_concentration"]

    # Separate properties and associated temperatures from data
    if prop_name == "Density":
        T = gcm.C2K(15)
    elif prop_name == "Viscosity":
        T = gcm.C2K(-20)
    
    # Vector for FuelLib predictions
    prop_pred = np.zeros_like(blend)

    if prop_name == "Density":
        for i in range(0, len(prop_pred)):
            # initial liquid mass fractions
            Y_li = blend[i]*fuel.Y_0 + (1-blend[i])*jetA.Y_0
            # Mixture density (returns rho in kg/m^3)
            prop_pred[i] = fuel.mixture_density(Y_li, T)
            # Convert density to CGS (g/cm^3)
            prop_pred[i] *= 1.0e-03

    if prop_name == "Viscosity":
        for i in range(0, len(prop_pred)):
            # initial liquid mass fractions
            Y_li = blend[i]*fuel.Y_0 + (1-blend[i])*jetA.Y_0
            # Mass of the droplet at current temp
            mass = gcm.drop_mass(fuel, drop["r_0"], Y_li, T)
            prop_pred[i] = fuel.mixture_kinematic_viscosity(mass, T)
            # Convert viscosity to mm^2/s
            prop_pred[i] *= 1.0e06

    return T, prop_data, blend_data, prop_pred


figW = 5.25 * len(prop_names)
fig, ax = plt.subplots(1, len(prop_names), figsize=(figW, 5.5), constrained_layout=True)

for i in range(len(prop_names)):

    for fuel_name in fuel_names:
        T, prop_data, blend_data, pred = getPredAndData(fuel_name, prop_names[i], blends)
        line_color, marker_style = linespecs(fuel_name)

        # Plot GCM predictions and data
        ax[i].plot(
            blends,
            pred,
            "-",
            color=line_color,
            label=leglab(fuel_name),
            linewidth=line_thickness,
        )

        ax[i].scatter(
            blend_data[1:],
            prop_data[1:],
            marker=marker_style,
            label=None,
            facecolors=line_color,
            s=marker_size,
        )

    # Data for pure Jet-A
    ax[i].scatter(
            blend_data[0],
            prop_data[0],
            marker="o",
            label="Jet-A",
            facecolors="darkorange",
            s=marker_size+2,
    )
    
    # Add labels and adjust ticks
    ax[i].set_xlabel("HEFA Concentration [wt %]", fontsize=fsize)
    ax[i].set_xticks([0, 20, 40, 60, 80, 100])
    ax[i].set_ylabel(ylab(prop_names[i],gcm.K2C(T)), fontsize=fsize)
    ax[i].tick_params(labelsize=ticksize)

handles, labels = ax[0].get_legend_handles_labels()
ax[i].legend(
    handles, labels, ncol=1, fontsize=fsize-2
)

plt.show()
