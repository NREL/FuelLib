import os
import sys
import numpy as np
import pandas as pd
import argparse
import FuelLib as fl

# Add the FuelLib directory to the Python path
FUELLIB_DIR = os.path.dirname(os.path.dirname(__file__))
if FUELLIB_DIR not in sys.path:
    sys.path.append(FUELLIB_DIR)
from paths import *

"""
Script that exports mixture properties over large temperature range for use in
Converge simulations.

This script is designed to be run from the command line and will create
a file named "mixturePropsGCM_<fuel_name>.csv" in the specified directory.
The file contains mixture properties for the fuel, formatted for Converge.

Usage:
    python Export4Converge.py --fuel_name <fuel_name>

Options:
        --units <units>
        --temp_min <temp_min> (K) 
        --temp_max <temp_max> (K)
        --temp_step <temp_step> (K)
        --export_dir <export_dir>
        --export_mix <export_mix>
"""


def export_converge(
    fuel,
    path=os.path.join(FUELLIB_DIR, "exportData"),
    units="mks",
    temp_min=0,
    temp_max=1000,
    temp_step=10,
    export_mix=False,
):
    """
    Export fuel properties to .csv for Converge simulations.

    :param fuel: An instance of the fuel class.
    :type fuel: fuel object

    :param path: Directory to save the input file.
    :type path: str, optional (default: FuelLib/exportData)

    :param units: Units for the properties ("mks" for SI, "cgs" for CGS).
    :type units: str, optional (default: "mks")

    :param temp_min: Minimum temperature (K) for the property calculations.
    :type temp_min: float, optional (default: 0)

    :param temp_max: Maximum temperature (K)for the property calculations.
    :type temp_max: float, optional (default: 1000)

    :param temp_step: Step size for temperature (K).
    :type temp_step: int, optional (default: 10)

    :export_mix: Whether to export individual component or mixture properties
    :type export_mix: bool, optional (default: False)

    :return: None
    :rtype: None
    """

    if not os.path.exists(path):
        os.makedirs(path)

    if export_mix:
        # Export mixture properties only
        file_name = os.path.join(path, f"mixturePropsGCM_{fuel.name}.csv")
        components = [fuel.name]
    else:
        # Export individual component properties and composition
        path = os.path.join(path, fuel.name)
        components = fuel.compounds

    # Unit conversion factors:
    if units.lower() == "cgs":
        # Convert from MKS to CGS
        conv_mw = 1e3  # kg/mol to g/mol
        conv_mu = 1e2  # Pa*s to Poise
        conv_surfacetension = 1e7  # N/m to dyne/cm
        conv_Lv = 1e4  # J/kg to erg/g
        conv_P = 1e1  # Pa to dyne/cm^2
        conv_rho = 1e3  # kg/m^3 to g/cm^3
        conv_Cl = 1e4  # J/kg/K to erg/g/K
        conv_thermcond = 1e5  # W/m/K to erg/cm/s/K
    else:
        conv_mw = 1
        conv_mu = 1
        conv_surfacetension = 1
        conv_Lv = 1
        conv_P = 1
        conv_rho = 1
        conv_Cl = 1
        conv_thermcond = 1

    # Assume droplet of 50 microns to account for compositional changes with temp
    drop_r = 50 * 1e-6  # initial droplet radius (m)

    # Vector of evenly space temperatures
    nT = int((temp_max - temp_min) / temp_step) + 1
    T = np.linspace(temp_min, temp_max, nT)

    # Round to nearest multiple of temp_step
    def nearest_temp(x, base=temp_step):
        return base * round(x / base)

    def nearest_floor(array, value):
        """
        Find the largest value in the array that is less than or equal to the given value.
        """
        if np.any(array <= value):
            return array[array <= value].max()
        else:
            raise ValueError(
                f"No temperature in the array is less than or equal to the critical point {value}. Choose a lower temp_min"
            )

    def nearest_ceil(array, value):
        """
        Find the smallest value in the array that is greater than or equal to the given value.
        """
        if np.any(array >= value):
            return array[array >= value].min()
        else:
            # Report an error if no value is found
            raise ValueError(
                f"No temperature in the array is greater than or equal the freezing point {value}. Choose a higher temp_max"
            )

    # Estimate freezing point and critical temp of mixture
    T_freeze = fl.mixing_rule(fuel.Tm, fuel.Y2X(fuel.Y_0))
    T_crit = fl.mixing_rule(fuel.Tc, fuel.Y2X(fuel.Y_0))
    T_min_allowed = nearest_temp(T_freeze)
    T_max_allowed = min(fuel.Tc)

    print(f"\nEstimated mixture freezing temp: {T_freeze:.2f} K")
    print(f"Min freezing temp min(Tm_i): {min(fuel.Tm):.2f} K")
    print(f"Max freezing temp max(Tm_i): {max(fuel.Tm):.2f} K")

    if np.any(T < T_min_allowed):
        T_min_allowed = nearest_ceil(T, T_min_allowed)
        # Set T_min_allowed to be the next temperature above T_min_allowed in T
        print(
            "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
        )
        print(
            f"   Warning: Some compounds have freezing temperatures above the estimated\n"
            f"   freezing temperature of the mixture ({T_freeze:.2f} K). All properties calculated\n"
            f"   below {T_min_allowed} will be set using a temperature of {T_min_allowed} K."
        )
        print(
            "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
        )

        print(f"\nEstimated mixture critical temp: {T_crit:.2f} K")
        print(f"Min critical temp min(Tc_i): {min(fuel.Tc):.2f} K")
        print(f"Max critical temp max(Tc_i): {max(fuel.Tc):.2f} K")
        if np.any(T > T_max_allowed):
            T_max_allowed = nearest_floor(T, T_max_allowed)
            print(
                "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
            )
            print(
                f"   Warning: Some compounds have critical temperatures below the estimated\n"
                f"   critical temperature of the mixture ({T_crit:.2f} K). All properties calculated\n"
                f"   above {T_max_allowed} will be set using a temperature of {T_max_allowed} K."
            )
            print(
                "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
            )

    for n, compound in enumerate(components):
        mu = np.zeros_like(T)  # Dynamic viscosity
        surface_tension = np.zeros_like(T)  # Surface tension
        Lv = np.zeros_like(T)  # Latent heat of vaporization
        pv = np.zeros_like(T)  # vapor pressure
        rho = np.zeros_like(T)  # density
        Cl = np.zeros_like(T)  # Specific heat
        thermal_conductivity = np.zeros_like(T)  # Thermal conductivity

        # Calculate GCM properties for a range of temperatures
        print(
            f"\nCalculating properties over {len(T)} temperatures from {temp_min} K to {temp_max} K..."
        )
        for k in range(len(T)):
            if T[k] <= T_min_allowed:
                Temp = T_min_allowed
            elif T[k] >= T_max_allowed:
                Temp = T_max_allowed
            else:
                Temp = T[k]

            if export_mix:
                Y_li = fuel.Y_0
                X_li = fuel.Y2X(Y_li)

                # Standard mixing rules for properties
                rho[k] = fuel.mixture_density(Y_li, Temp)  # kg/m^3
                mu[k] = fuel.mixture_dynamic_viscosity(Y_li, Temp)  # Pa*s
                pv[k] = fuel.mixture_vapor_pressure(Y_li, Temp)  # Pa
                surface_tension[k] = fuel.mixture_surface_tension(Y_li, Temp)  # N/m
                thermal_conductivity[k] = fuel.mixture_thermal_conductivity(Y_li, Temp)

                # Generic mixing rules for latent heat and specific heat
                Lv[k] = fl.mixing_rule(
                    fuel.latent_heat_vaporization(Temp), X_li
                )  # J/kg
                Cl[k] = fl.mixing_rule(fuel.Cl(Temp), X_li)  # J/kg/K
            else:
                rho[k] = fuel.density(Temp)[n]  # kg/m^3
                mu[k] = fuel.viscosity_dynamic(Temp)[n]  # Pa*s
                pv[k] = fuel.psat(Temp)[n]  # Pa
                surface_tension[k] = fuel.surface_tension(Temp)[n]  # N/m
                thermal_conductivity[k] = fuel.thermal_conductivity(Temp)[n]
                Lv[k] = fuel.latent_heat_vaporization(Temp)[n]  # J/kg
                Cl[k] = fuel.Cl(Temp)[n]  # J/kg/K

        if units.lower() == "cgs":
            # Convert properties to CGS units
            data = {
                "Temperature (K)": T,
                "Critical Temperature (K)": T_crit + np.zeros_like(T),
                "Viscosity (Poise)": mu * conv_mu,
                "Surface Tension (dyne/cm)": surface_tension * conv_surfacetension,
                "Heat of Vaporization (erg/g)": Lv * conv_Lv,
                "Vapor Pressure (dyne/cm^2)": pv * conv_P,
                "Density (g/cm^3)": rho * conv_rho,
                "Specific Heat (erg/g/K)": Cl * conv_Cl,
                "Thermal Conductivity (erg/cm/s/K)": thermal_conductivity
                * conv_thermcond,
            }

        else:
            # MKS units
            data = {
                "Temperature (K)": T,
                "Critical Temperature (K)": T_crit + np.zeros_like(T),
                "Viscosity (Pa*s)": mu,
                "Surface Tension (N/m)": surface_tension,
                "Heat of Vaporization (J/kg)": Lv,
                "Vapor Pressure (Pa)": pv,
                "Density (kg/m^3)": rho,
                "Specific Heat (J/kg/K)": Cl,
                "Thermal Conductivity (W/m/K)": thermal_conductivity,
            }

        # Write the properties to the input file
        if export_mix:
            print(f"\nWriting mixture properties to {file_name}")
        else:
            file_name = os.path.join(path, f"{n}_{compound}.csv")
            print(f"\nWriting properties for {compound} to {file_name}")
        if not os.path.exists(path):
            os.makedirs(path)
        if os.path.exists(file_name):
            os.remove(file_name)
        df = pd.DataFrame(data)
        df.to_csv(file_name, index=False)

    if not export_mix:
        # Also export the initial mass fractions
        file_name = os.path.join(path, f"composition_{fuel.name}.csv")
        if os.path.exists(file_name):
            os.remove(file_name)
        print(f"\nWriting mass fractions for {fuel.name} to {file_name}")
        if units.lower() == "cgs":
            mw_label = "Molecular Weight (g/mol)"
        else:
            mw_label = "Molecular Weight (kg/mol)"
        data = {
            "Index": range(len(fuel.compounds)),
            "Component": fuel.compounds,
            "Mass Fraction": fuel.Y_0,
            "Mole Fraction": fuel.Y2X(fuel.Y_0),
            mw_label: fuel.MW * conv_mw,
        }
        dfY = pd.DataFrame(data)
        dfY.to_csv(file_name, index=False)


def main():
    """
    Main function to execute the export process.

    :param --fuel_name: Name of the fuel (mandatory).
    :type --fuel_name: str

    :param --fuel_data_dir: Directory where fuel data files are located.
    :type --fuel_data_dir: str, optional (default: FuelLib/fuelData)

    :param --units: Units for critical properties. Options are mks or cgs.
    :type --units: str, optional (default: mks)

    :param --temp_min: Minimum temperature (K) for the property calculations.
    :type --temp_min: float, optional (default: 0 K)

    :param --temp_max: Maximum temperature (K) for the property calculations.
    :type --temp_max: float, optional (default: 1000 K)

    :param --temp_step: Step size for temperature (K).
    :type --temp_step: float, optional (default: 10 K)

    :param --export_dir: Directory to export the properties.
    :type --export_dir: str, optional (default: FuelLib/exportData)

    :param --export_mix: Whether to export individual component or mixture properties.
    :type --export_mix: bool, optional (default: False)

    :raises FileNotFoundError: If required files for the specified fuel are not found.
    """

    # Set up argument parser
    parser = argparse.ArgumentParser(
        description="Export mixture fuel properties for Converge simulations."
    )

    # Mandatory argument for fuel name
    parser.add_argument(
        "--fuel_name",
        required=True,
        help="Name of the fuel (mandatory).",
    )

    # Optional argument for fuel data directory
    parser.add_argument(
        "--fuel_data_dir",
        default=FUELDATA_DIR,
        help="Directory where fuel data files are located (optional, default: FuelLib/fuelData).",
    )

    # Optional argument for units
    # Default is 'mks', but can be set to 'cgs'
    parser.add_argument(
        "--units",
        default="mks",
        help="Units for critical properties: mks or cgs (optional, default: mks).",
    )

    # Optional argument for minimum temperature
    parser.add_argument(
        "--temp_min",
        type=float,
        default=0,
        help="Minimum temperature (K) for the property calculations (optional, default: 0).",
    )

    # Optional argument for maximum temperature
    parser.add_argument(
        "--temp_max",
        type=float,
        default=1000,
        help="Maximum temperature (K) for the property calculations (optional, default: 1000).",
    )

    # Optional argument for temperature step size
    parser.add_argument(
        "--temp_step",
        type=int,
        default=10,
        help="Step size for temperature (K) (optional, default: 10).",
    )

    # Optional argument for export directory
    parser.add_argument(
        "--export_dir",
        default=os.path.join(FUELLIB_DIR, "exportData"),
        help="Directory to export the properties (optional, default: FuelLib/exportData).",
    )

    # Optional argument for exporting mixture properties
    parser.add_argument(
        "--export_mix",
        type=lambda x: str(x).lower() in ["true", "1"],
        default=False,
        help="Option to export mixture properties of the fuel (True or False, default: False).",
    )

    # Parse arguments
    args = parser.parse_args()
    fuel_name = args.fuel_name
    fuel_data_dir = args.fuel_data_dir
    units = args.units.lower()
    temp_min = args.temp_min
    temp_max = args.temp_max
    temp_step = args.temp_step
    export_dir = args.export_dir
    export_mix = args.export_mix

    # Print the parsed arguments
    print(f"Preparing to export mixture properties:")
    print(f"    Fuel name: {fuel_name}")
    if export_mix:
        print(f"    Exporting mixture properties: True")
    print(f"    Units: {units}")
    print(f"    Minimum temperature: {temp_min} K")
    print(f"    Maximum temperature: {temp_max} K")
    print(f"    Temperature step size: {temp_step} K")
    print(f"    Export directory: {export_dir}")
    print(f"    Fuel data directory: {fuel_data_dir}")

    # Check if necessary files exist in the fuelData directory
    print("\nChecking for required files...")
    gcxgc_file = os.path.join(fuel_data_dir, f"gcData/{fuel_name}_init.csv")
    decomp_file = os.path.join(fuel_data_dir, f"groupDecompositionData/{fuel_name}.csv")
    if not os.path.exists(gcxgc_file):
        err = f"GCXGC file for {fuel_name} not found in {fuel_data_dir}/gcData. gcxgc_file = {gcxgc_file}"
        raise FileNotFoundError(err)
    if not os.path.exists(decomp_file):
        err = f"Decomposition file for {fuel_name} not found in {fuel_data_dir}/groupDecompositionData. decomp_file = {decomp_file}"
        raise FileNotFoundError(err)
    print("All required files found.")

    # Create the groupContribution object for the specified fuel
    fuel = fl.fuel(fuel_name, fuelDataDir=fuel_data_dir)

    # Export properties for Pele
    export_converge(
        fuel,
        path=export_dir,
        units=units,
        temp_min=temp_min,
        temp_max=temp_max,
        temp_step=temp_step,
        export_mix=export_mix,
    )

    print("\nExport completed successfully!")


if __name__ == "__main__":
    main()
