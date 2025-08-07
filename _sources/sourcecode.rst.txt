Source Code
===========

This page provides an overview of the source code available at `github.com/NREL/FuelLib <https://github.com/NREL/FuelLib>`_.

.. _source-code-structure:

FuelLib File Organization
-------------------------

- **docs:** directory containing the documentation source files
- **examples:** directory containing example scripts that demonstrate how to use FuelLib

    - *ex_compositionPlots.py*: example script that generates composition plots for a given fuel
    - *ex_hefaBlends.py*: example script that calculates properties of HEFA:Jet-A blends
    - *ex_mixtureProperties.py*: validation script that calculates properties of single component fuels and mixture properties of multicomponent fuels.
- **fuelData:** 
    - **gcData:** directory containing a collection of GCxGC compositional data by weight percentages
    - **groupDecompositionData:** directory containing a collection of functional group decompositions
    - **propertiesData:** directory containing measurement or predicted data for validation (see *fuelData/dataReferences.md*)
- **gcmTableData:** directory that contains the pre-tabulated group contributions
- *GroupContributionMethod.py*: class for enabling GCM predictions
- **tests:**  directory containing CI unit tests for FuelLib. The CI test checks if the cumulative error of property predictions of a new proposed model are less than or equal to the current model.
    
    - **baselinePredictions:** directory that contains baseline predictions
    - *test_accuracy.py*: unit test used in CI for verifying new model predictions preserve accuracy
    - *test_baseline.py*: generates .csv files for the baseline model predictions, which are stored in **baselinePredictions**
    - *test_functions.py*: collection of functions used by *test_baseline.py* and *test_accuracy.py*.   


Source Code Auto-Documentation
------------------------------
Click on GroupContributionMethod below for the full auto-documentation.

.. autosummary::
    :toctree: generated

    GroupContributionMethod