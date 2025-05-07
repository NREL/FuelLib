import os
import numpy as np
import pandas as pd
import test_functions as fxns
import unittest


class CompTestCase(unittest.TestCase):
    """Test that accuracy is preserved"""

    def test_accuracy(self):
        """Does the PR impact the accuracy of single component fuel predictions?"""
        fuellib_dir = os.path.dirname(__file__)
        baseline_dir = os.path.join(fuellib_dir, "fuelData/baselinePredictions")

        # Set the precentage for variation in max error
        max_error_diff = 1e-6

        # Fuels to test
        fuel_names = [
            "heptane",
            "decane",
            "dodecane",
            "posf10264",
            "posf10325",
            "posf10289",
        ]

        # Properties to test
        prop_names = [
            "Density",
            "Viscosity",
            "VaporPressure",
            "SurfaceTension",
            "ThermalConductivity",
        ]

        # droplet specs
        drop = {}
        drop["d_0"] = (
            100 * 1e-6
        )  # initial droplet diameter (m), note: size doesn't matter
        drop["r_0"] = drop["d_0"] / 2.0  # initial droplet radius (m)
        drop["V_0"] = 4.0 / 3.0 * np.pi * drop["r_0"] ** 3  # initial droplet volume

        # Compare to NIST predictions and previous model predictions
        for fuel_name in fuel_names:

            baseline_file = os.path.join(baseline_dir, f"{fuel_name}.csv")
            df_base = pd.read_csv(baseline_file, skiprows=[1])

            for prop in prop_names:
                # Baseline error
                err_base = df_base[f"Error_{prop}"].dropna().to_numpy()
                # err_base += err_base*max_error_diff

                # Cumulative baseline error
                sum_err_base = np.sum(err_base)

                # Allow cumulative error within max_error_diff
                sum_err_base += sum_err_base * max_error_diff

                # Get predictions for current model
                T, data, pred = fxns.getPredAndData(drop, fuel_name, prop)
                err = np.abs(data - pred)
                sum_err = np.sum(err)

                np.testing.assert_array_less(
                    sum_err,
                    sum_err_base,
                    err_msg=f"Prediction of {prop} for {fuel_name} is less accurate than previous model",
                )


if __name__ == "__main__":
    unittest.main()
