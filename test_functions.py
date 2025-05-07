import os
import numpy as np
import pandas as pd
import GroupContributionMethod as gcm

def drop_mass(drop, fuel, Yi, T):
    return drop["V_0"] / (fuel.molar_liquid_vol(T) @ Yi) * Yi * fuel.MW  # (kg)


def getPredAndData(drop, fuel_name, prop_name):
    # Get the fuel properties based on the GCM
    fuel = gcm.groupContribution(fuel_name)

    # initial liquid mass fractions
    Y_li = fuel.Y_0

    data_file = f"{fuel_name}.csv"
    dataPath = os.path.join(fuel.fuelDataDir, "propertiesData")
    data = pd.read_csv(os.path.join(dataPath, data_file), skiprows=[1])

    # Separate properties and associated temperatures from data
    T_data = data.Temperature[data[prop_name].notna()].to_numpy(dtype=float)
    prop_data = data[prop_name].dropna().to_numpy()

    # Vector for predictions
    pred = np.zeros_like(T_data)

    T_pred = gcm.C2K(T_data)

    if prop_name == "Density":
        for i in range(0, len(T_pred)):
            # Mixture density (returns rho in kg/m^3)
            pred[i] = fuel.mixture_density(fuel.Y_0, T_pred[i])
            # Convert density to CGS (g/cm^3)
            pred[i] *= 1.0e-03

    if prop_name == "VaporPressure":
        for i in range(0, len(T_pred)):
            # Mass of the droplet at current temp
            mass = drop_mass(drop, fuel, Y_li, T_pred[i])
            # Mixture vapor pressure (returns pv in Pa)
            pred[i] = fuel.mixture_vapor_pressure(mass, T_pred[i])
            # Convert vapor pressure to kPa
            pred[i] *= 1.0e-03

    if prop_name == "Viscosity":
        for i in range(0, len(T_pred)):
            # Mass of the droplet at current temp
            mass = drop_mass(drop, fuel, Y_li, T_pred[i])
            pred[i] = fuel.mixture_kinematic_viscosity(mass, T_pred[i])
            # Convert viscosity to mm^2/s
            pred[i] *= 1.0e06

    if prop_name == "SurfaceTension":
        for i in range(0, len(T_pred)):
            # Mass of the droplet at current temp
            mass = drop_mass(drop, fuel, Y_li, T_pred[i])
            pred[i] = fuel.mixture_surface_tension(mass, T_pred[i])

    if prop_name == "ThermalConductivity":
        for i in range(0, len(T_pred)):
            # Mass of the droplet at current temp
            mass = drop_mass(drop, fuel, Y_li, T_pred[i])
            pred[i] = fuel.mixture_thermal_conductivity(mass, T_pred[i])

    return T_data, prop_data, pred