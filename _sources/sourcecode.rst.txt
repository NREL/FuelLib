Source Code
===========

This page provides an overview of the source code available at `github.com/NREL/FuelLib <https://github.com/NREL/FuelLib>`_.

.. _source-code-structure:

FuelLib File Structure
----------------------

- *GroupContributionMethod.py*: class for enabling GCM predictions
- *ex_mixtureProperties.py*: validation script that calculates mixture properties
- *ex_hefaBlends.py*: example script that calculates properties of HEFA:Jet-A blends
- **gcmTableData:** directory that contains the pre-tabulated group contributions
- **fuelData:** 
    - **gcData:** directory containing a collection of GCxGC compositional data by weight percentages
    - **groupDecompositionData:** directory containing a collection of functional group decompositions
    - **propertiesData:** directory containing measurement or predicted data for validation (see *dataReferences.md*)
    - **baselinePredictions:** directory that contains FuelLib's baseline predictions for specified fuels and properties, which are used in CI test to ensure future changes preserve model accuracy (see next section).

CI testing scripts
^^^^^^^^^^^^^^^^^^
The following scripts are used to ensure new changes don't negatively impact FuelLib's
predictive capabilities.  Specifically, the CI test checks if the cumulative error of 
property predictions of a new proposed model are less than or equal to the current model.

- *test_baseline.py*: generates .csv files for the baseline model predictions, which are stored in **baselinePredictions**
- *test_accuracy.py*: unit test for verifying new model predictions preserve accuracy
- *test_functions.py*: collection of functions used by *test_baseline.py* and *test_accuracy.py*.   

API
---
Click on GroupContributionMethod below for the full API. 

.. autosummary::
    :toctree: generated

    GroupContributionMethod