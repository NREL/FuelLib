import os
import sys

FUELLIB_DIR = os.path.dirname(__file__)
GCMTABLE_DIR = os.path.join(FUELLIB_DIR, "gcmTableData")
SOURCE_DIR = os.path.join(FUELLIB_DIR, "source")
FUELDATA_DIR = os.path.join(FUELLIB_DIR, "fuelData")
FUELDATA_GC_DIR = os.path.join(FUELDATA_DIR, "gcData")
FUELDATA_DECOMP_DIR = os.path.join(FUELDATA_DIR, "groupDecompositionData")
FUELDATA_PROPS_DIR = os.path.join(FUELDATA_DIR, "propertiesData")
TESTS_DIR = os.path.join(FUELLIB_DIR, "tests")
TESTS_BASELINE_DIR = os.path.join(TESTS_DIR, "baselinePredictions")
TUTORIALS_DIR = os.path.join(FUELLIB_DIR, "tutorials")

sys.path.append(SOURCE_DIR)
