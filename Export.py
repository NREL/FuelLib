import pandas as pd
import os
import sys
import GroupContributionMethod as gcm
"""
Script that exports critical properties and initial mass and mole fraction data
for use in Pele simulations.

Usage: 
python Export.py [fuel_name]
"""

def export_4_pele(fuel, path="FuelLibProps", units="mks"):
    """
    Export fuel properties to .cpp and .H files for use in Pele.

    :param fuel: An instance of the groupContribution class.
    :type fuel: groupContribution object
    :param path: Directory to save the CSV file.
    :type path: str, optional
    :param units: Units for the properties ("mks" for SI, "cgs" for CGS).
    :type units: str, optional
    :return: None
    :rtype: None
    """
    
    if not os.path.exists(path):
        os.makedirs(path)
    

    # Terms for liquid specific heat capacity in (J/kg/K) or (erg/g/K)
    Cp_stp = fuel.Cp_stp / fuel.MW  
    Cp_B = fuel.Cp_B / fuel.MW  
    Cp_C = fuel.Cp_C / fuel.MW  

    # Conversion factors:
    if units.lower() == "cgs":
        # Convert from MKS to CGS
        conv_MW = 1e3 # kg/mol to g/mol
        conv_Cp = 1e7 # J/kg/K to erg/g/K
        conv_Vm = 1e6 # m^3/mol to cm^3/mol
        conv_Lv = 1e3 # J/kg to erg/g
        conv_P  = 1e1 # Pa to dyne/cm^2
    else:
        conv_MW = 1e3 # kg/mol to g/mol
        conv_Cp = 1.0 
        conv_Vm = 1.0 
        conv_Lv = 1.0 
        conv_P  = 1.0

    df = pd.DataFrame({
        "Compound": fuel.name,
        "MW": fuel.MW * conv_MW, 
        "Tc": fuel.Tc, 
        "Pc": fuel.Pc * conv_P,
        "Vc": fuel.Vc,
        "Tb": fuel.Tb,
        "omega": fuel.omega,
        "Vm_stp": fuel.Vm_stp,
        "Cp_stp": Cp_stp,
        "Cp_B": Cp_B,
        "Cp_C": Cp_C,
        "Lv_stp": fuel.Lv_stp
    })
    
    # Create .H and .cpp files for export
    header_file = os.path.join(path, "FuelLibProps.H")
    cpp_file = os.path.join(path, "FuelLibProps.cpp")
    with open(header_file, 'w') as f:
        f.write("#ifndef FUELLIBPROPS_H\n#define FUELLIBPROPS_H\n")
        f.write("\n#include <AMReX_Gpu.H>\n#include <AMReX_REAL.H>\n")
    with open(cpp_file, 'w') as f:
        f.write("#include FuelLibProps.H\n")
    
    # Write compound names
    
    # Close files
    with open(header_file, 'a') as f:
        f.write("#endif // FUELLIBPROPS_H\n")

def main():
    """
    Main function to execute the export process.
    """

    # For testing and debugging purposes, use a default fuel name
    fuel_name = "posf10325"

    # Export directory given by second command line argument or default to "FuelLibProps"
    export_dir = "FuelLibProps"
    #if len(sys.argv) > 1:
    #    fuel_name = sys.argv[1]
    #    if len(sys.argv) > 2:
    #        export_dir = sys.argv[2]

    #if len(sys.argv) > 1:
    #fuel_name = sys.argv[1]
    # Check if necessary files exist in the fuelData directory
    decomp_dir = os.path.join(gcm.groupContribution.fuelDataDir, "groupDecompositionData")
    gcxgc_dir = os.path.join(gcm.groupContribution.fuelDataDir, "gcData")
    gcxgc_file = os.path.join(gcxgc_dir, f"{fuel_name}_init.csv")
    decomp_file = os.path.join(decomp_dir, f"{fuel_name}.csv")
    if not os.path.exists(gcxgc_file):
        raise FileNotFoundError(f"GCXGC file for {fuel_name} not found in {gcxgc_dir}. gxcgc_file = {gcxgc_file}")
    if not os.path.exists(decomp_file):
        raise FileNotFoundError(f"Decomposition file for {fuel_name} not found in {decomp_dir}.")
    #else: 
    #    # Throw an error if no fuel name is provided
    #    raise ValueError("Please provide a fuel name as a command line argument.")

    fuel = gcm.groupContribution(fuel_name)

    # Export properties for Pele
    export_4_pele(fuel, path=export_dir, units="mks")   

if __name__ == "__main__":
    main()