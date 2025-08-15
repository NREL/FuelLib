import os
import numpy as np
import pandas as pd
import sys

# Add the FuelLib directory to the Python path
fuellib_dir = os.path.dirname(os.path.dirname(__file__))
sys.path.append(fuellib_dir)
import GroupContributionMethod as gcm


def getPredAndData(fuel_name, prop_name):
    # Get the fuel properties based on the GCM
    fuel = gcm.groupContribution(fuel_name)

    # Assume a droplet of fuel with radius 50 microns (r doesn't matter)
    drop_r = 50 * 1e-6  # initial droplet diameter (m)

    data_file = f"{fuel_name}.csv"
    dataPath = os.path.join(fuel.fuelDataDir, "propertiesData")
    data = pd.read_csv(os.path.join(dataPath, data_file), skiprows=[1])

    # Separate properties and associated temperatures from data
    T_data = data.Temperature[data[prop_name].notna()].to_numpy(dtype=float)
    prop_data = data[prop_name].dropna().to_numpy()

    # Vector for predictions
    pred = np.zeros_like(T_data)

    # Vectors for temperature (convert from C to K)
    T_pred = gcm.C2K(T_data)

    for i in range(0, len(T_pred)):
        # Correct droplet mass (GCxGC at standard temperature)
        mass = gcm.droplet_mass(fuel, drop_r, fuel.Y_0, T_pred[i])
        Y_li = fuel.mass2Y(mass)

        if prop_name == "Density":
            # Mixture density (returns rho in kg/m^3)
            pred[i] = fuel.mixture_density(Y_li, T_pred[i])
            # Convert density to CGS (g/cm^3)
            pred[i] *= 1.0e-03

        if prop_name == "VaporPressure":
            # Mixture vapor pressure (returns pv in Pa)
            pred[i] = fuel.mixture_vapor_pressure(Y_li, T_pred[i])
            # Convert vapor pressure to kPa
            pred[i] *= 1.0e-03

        if prop_name == "Viscosity":
            pred[i] = fuel.mixture_kinematic_viscosity(Y_li, T_pred[i])
            # Convert viscosity to mm^2/s
            pred[i] *= 1.0e06

        if prop_name == "SurfaceTension":
            pred[i] = fuel.mixture_surface_tension(Y_li, T_pred[i])

        if prop_name == "ThermalConductivity":
            pred[i] = fuel.mixture_thermal_conductivity(Y_li, T_pred[i])

    return T_data, prop_data, pred
