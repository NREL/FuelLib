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
Script that exports critical properties and initial mass fraction data
for use in Pele simulations.

This script is designed to be run from the command line and will create
a file named "sprayPropsGCM_<fuel_name>.inp" in the specified directory.
The file contains properties for each compound in the fuel, formatted for Pele.

Usage:
    python Export4Pele.py --fuel_name <fuel_name>

Options:
        --units <units>
        --dep_fuel_names <dep_fuel_names>
        --max_dep_fuels <max_dep_fuels>
        --export_dir <export_dir>
"""


def vec_to_str(vec):
    """
    Convert a list or numpy array to a string representation.

    :param vec: List or numpy array to convert.
    :return: String representation of the vector.
    """

    # If strings return string[0] string[1] ... string[n]
    if isinstance(vec, list):
        return " ".join(f"{v}" for v in vec)
    # Else if numbers, format with spaces between no commas or []
    elif isinstance(vec, (pd.Series, pd.DataFrame)):
        return " ".join(f"{v}" for v in vec.values)


def export_pele(
    fuel,
    path=os.path.join(FUELLIB_DIR, "exportData"),
    units="mks",
    dep_fuel_names=None,
    max_dep_fuels=30,
):
    """
    Export fuel properties to input file for Pele simulations.

    :param fuel: An instance of the fuel class.
    :type fuel: fuel object

    :param path: Directory to save the input file.
    :type path: str, optional (default: FuelLib/exportData)

    :param units: Units for the properties ("mks" for SI, "cgs" for CGS).
    :type units: str, optional (default: "mks")

    :param dep_fuel_names: List or single fuel that each compound deposits to.
    :type dep_fuel_names: list of str, optional (default: None)

    :param max_dep_fuels: Maximum number of deposition fuels to consider.
    :type max_dep_fuels: int, optional (default: 30)

    :return: None
    :rtype: None
    """

    if not os.path.exists(path):
        os.makedirs(path)

    # Names of the input file
    file_name = os.path.join(path, f"sprayPropsGCM_{fuel.name}.inp")

    # If dep_fuel_names is not provided, use fuel.compounds
    if dep_fuel_names is None:
        if len(fuel.compounds) <= max_dep_fuels:
            # If no deposition fuel names are provided, use the compounds as deposition fuels
            dep_fuel_names = fuel.compounds
        else:
            # If more than max_dep_fuels, deposit all compoudns to fuel.name.upper()
            # This assumes a POSF fuel with a single deposition fuel
            dep_fuel_names = [fuel.name.upper()] * len(fuel.compounds)
    elif len(dep_fuel_names) == 1:
        # If a single deposition fuel name is provided, use it for all compounds
        dep_fuel_names = [dep_fuel_names[0]] * len(fuel.compounds)
    elif len(dep_fuel_names) != len(fuel.compounds):
        raise ValueError(
            "Length of dep_fuel_names must be one or match the number of compounds in the fuel."
        )

    # Unit conversion factors:
    if units.lower() == "cgs":
        # Convert from MKS to CGS
        conv_MW = 1e3  # kg/mol to g/mol
        conv_Cp = 1e4  # J/kg/K to erg/g/K
        conv_Vm = 1e6  # m^3/mol to cm^3/mol
        conv_Lv = 1e4  # J/kg to erg/g
        conv_P = 1e1  # Pa to dyne/cm^2
    else:
        conv_MW = 1.0
        conv_Cp = 1.0
        conv_Vm = 1.0
        conv_Lv = 1.0
        conv_P = 1.0

    # Terms for liquid specific heat capacity in (J/kg/K) or (erg/g/K)
    # Cp(T) = Cp_A + Cp_B * theta + Cp_C * theta^2
    # where theta = (T - 298.15) / 700
    Cp_stp = fuel.Cp_stp / fuel.MW
    Cp_B = fuel.Cp_B / fuel.MW
    Cp_C = fuel.Cp_C / fuel.MW

    # Dataframe of all properties with unit conversions to be exported
    print("\nCalculating GCM properties at standard conditions...")
    df = pd.DataFrame(
        {
            "Compound": fuel.compounds,
            "Y_0": fuel.Y_0,
            "MW": fuel.MW * conv_MW,
            "Tc": fuel.Tc,
            "Pc": fuel.Pc * conv_P,
            "Vc": fuel.Vc * conv_Vm,
            "Tb": fuel.Tb,
            "omega": fuel.omega,
            "Vm_stp": fuel.Vm_stp * conv_Vm,
            "Cp_stp": Cp_stp * conv_Cp,
            "Cp_B": Cp_B * conv_Cp,
            "Cp_C": Cp_C * conv_Cp,
            "Lv_stp": fuel.Lv_stp * conv_Lv,
        }
    )
    # Get the property names
    prop_names = ["MW", "Tc", "Pc", "Vc", "Tb", "omega", "Vm_stp", "Cp_stp", "Cp_B", "Cp_C", "Lv_stp"]

    formatted_names = {
        "MW": ("molar_weight", ["kg/mol", "g/mol"]),
        "Tc": ("crit_temp", ["K", "K"]),
        "Pc": ("crit_press", ["Pa", "dyne/cm^2"]),
        "Vc": ("crit_vol", ["m^3/mol", "cm^3/mol"]),
        "Tb": ("boil_temp", ["K", "K"]),
        "omega": ("acentric_factor", ["-", "-"]),
        "Vm_stp": ("molar_vol", ["m^3/mol", "cm^3/mol"]),
        "Cp_stp": ("cp_a", ["J/kg/K", "erg/g/K"]),
        "Cp_B": ("cp_b", ["J/kg/K", "erg/g/K"]),
        "Cp_C": ("cp_c", ["J/kg/K", "erg/g/K"]),
        "Lv_stp": ("latent", ["J/kg", "erg/g"]),
    }

    # Get date and time for the header
    from datetime import datetime

    now = datetime.now()
    dt_string = now.strftime("%Y-%m-%d %H:%M:%S")

    # Get git commit hash for the header
    import subprocess

    try:
        git_commit = (
            subprocess.check_output(["git", "rev-parse", "HEAD"])
            .strip()
            .decode("utf-8")
        )
    except Exception:
        git_commit = "N/A"

    # Get the remote url for the header
    try:
        git_remote = (
            subprocess.check_output(["git", "config", "--get", "remote.origin.url"])
            .strip()
            .decode("utf-8")
        )
    except Exception:
        git_remote = "N/A"

    # Write the properties to the input file
    print(f"Writing properties to {file_name}...")
    if os.path.exists(file_name):
        os.remove(file_name)
    with open(file_name, "a") as f:
        f.write(
            f"# -----------------------------------------------------------------------------\n"
            f"# sprayPropsGCM_{fuel.name}.inp\n"
            f"# Generated on {dt_string}\n"
            f"# FuelLib remote URL: {git_remote}\n"
            f"# Git commit: {git_commit}\n"
            f"# Units: {units.upper()}\n"
            f"# -----------------------------------------------------------------------------\n\n"
        )
        f.write(f"particles.spray_fuel_num = {len(fuel.compounds)}\n")
        f.write(f"particles.fuel_species = {vec_to_str(df['Compound'].tolist())}\n")
        f.write(f"particles.Y_0 = {vec_to_str(df['Y_0'].tolist())}\n")
        f.write(f"particles.dep_fuel_names = {vec_to_str(dep_fuel_names)}\n")

        for comp_name in fuel.compounds:
            f.write(f"\n# Properties for {comp_name} in {units.upper()}\n")
            for prop in prop_names:
                if prop in formatted_names:
                    value = df.loc[df["Compound"] == comp_name, prop].values[0]
                    prop_name, unit_txt = formatted_names[prop]
                    if units.lower() == "cgs":
                        unit_txt = unit_txt[1]
                    else:
                        unit_txt = unit_txt[0]
                    # Write the property to the file
                    f.write(
                        f"particles.{comp_name}_{prop_name} = {value:.6f} # {unit_txt}\n"
                    )


def main():
    """
    Main function to execute the export process.

    :param --fuel_name: Name of the fuel (mandatory).
    :type --fuel_name: str

    :param --units: Units for critical properties. Options are "mks" (default) or "cgs".
    :type --units: str, optional

    :param --dep_fuel_names: Space-separated list with len(fuel.compounds) or single fuel that all compounds deposit. Default is fuel.compounds.
    :type --dep_fuel_names: str, optional

    :param --max_dep_fuels: Maximum number of deposition fuels to consider. Default is 30.
    :type --max_dep_fuels: int, optional

    :param --export_dir: Directory to export the properties. Default is "FuelLib/exportData".
    :type --export_dir: str, optional

    :raises FileNotFoundError: If required files for the specified fuel are not found.
    """

    # Set up argument parser
    parser = argparse.ArgumentParser(
        description="Export fuel properties for Pele simulations."
    )

    # Mandatory argument for fuel name
    parser.add_argument(
        "--fuel_name",
        required=True,
        help="Name of the fuel (mandatory).",
    )

    # Optional argument for units
    # Default is 'mks', but can be set to 'cgs'
    parser.add_argument(
        "--units",
        default="mks",
        help="Units for critical properties: mks or cgs (optional, default: mks).",
    )

    # Optional argument for deposition fuel names
    parser.add_argument(
        "--dep_fuel_names",
        nargs="+",  # Accepts one or more values
        default=None,
        help="Space-separated list or single fuel that each compound deposits to (optional, default: fuel.compounds).",
    )

    # Optional argument for maximum number of deposition fuels
    parser.add_argument(
        "--max_dep_fuels",
        type=int,
        default=30,
        help="Maximum number of deposition fuels to consider (optional, default: 30).",
    )

    # Optional argument for export directory
    parser.add_argument(
        "--export_dir",
        default=os.path.join(FUELLIB_DIR, "exportData"),
        help="Directory to export the properties (optional, default: FuelLib/exportData).",
    )

    # Parse arguments
    args = parser.parse_args()
    fuel_name = args.fuel_name
    units = args.units.lower()
    dep_fuel_names = args.dep_fuel_names
    max_dep_fuels = args.max_dep_fuels
    export_dir = args.export_dir

    # Print the parsed arguments
    print(f"Preparing to export properties:")
    print(f"    Fuel name: {fuel_name}")
    print(f"    Units: {units}")
    print(f"    Export directory: {export_dir}")

    # Check if necessary files exist in the fuelData directory
    print("\nChecking for required files...")
    gcxgc_file = os.path.join(FUELDATA_GC_DIR, f"{fuel_name}_init.csv")
    decomp_file = os.path.join(FUELDATA_DECOMP_DIR, f"{fuel_name}.csv")
    if not os.path.exists(gcxgc_file):
        raise FileNotFoundError(
            f"GCXGC file for {fuel_name} not found in {FUELDATA_GC_DIR}. gxcgc_file = {gcxgc_file}"
        )
    if not os.path.exists(decomp_file):
        raise FileNotFoundError(
            f"Decomposition file for {fuel_name} not found in {FUELDATA_DECOMP_DIR}."
        )
    print("All required files found.")

    # Create the groupContribution object for the specified fuel
    fuel = fl.fuel(fuel_name)

    # Export properties for Pele
    export_pele(
        fuel,
        path=export_dir,
        units=units,
        dep_fuel_names=dep_fuel_names,
        max_dep_fuels=max_dep_fuels,
    )

    print("\nExport completed successfully!")


if __name__ == "__main__":
    main()
