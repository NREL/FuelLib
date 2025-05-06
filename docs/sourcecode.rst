Source Code
===========

This page provides an overview of the source code available at `github.com/NREL/FuelLib <https://github.com/NREL/FuelLib>`_.

.. _source-code-structure:

FuelLib File Structure
----------------------

- *GroupContributionMethod.py*: class for enabling GCM predictions
- *ex-mixtureProperties.py*: validation script that calculates mixture properties
- **gcmTableData:** directory that contains the pre-tabulated group contributions
- **fuelData:** 
    - **gcData:** directory containing a collection of GCxGC compositional data by weight percentages
    - **groupDecompositionData:** directory containing a collection of functional group decompositions
    - **propertiesData:** directory containing measurement or predicted data for validation (see *dataReferences.md*)

API
---
Click on GroupContributionMethod below for the full API. 

.. autosummary::
    :toctree: generated

    GroupContributionMethod