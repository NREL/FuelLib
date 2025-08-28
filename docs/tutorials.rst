Tutorials
=========

This section provides tutorials on how to use the FuelLib library. 

Introduction
------------

Clone the FuelLib repository from GitHub: ::

    git clone https://github.com/NREL/FuelLib.git

Create and activate a Conda environment, install the required dependencies: ::

    conda create --name fuellib-env matplotlib pandas scipy
    conda activate fuellib-env

Change to the FuelLib/tutorials directory: ::

    cd FuelLib/tutorials

Required Input files
^^^^^^^^^^^^^^^^^^^^^

FuelLib requires two input files for any given fuel, ``<fuel_name>``:

- ``FuelLib/fuelData/gcData/<fuel_name>_init.csv``: an initial weight percentage composition of the fuel components
- ``FuelLib/fuelData/groupDecompositionData/<fuel_name>.csv``: the fundamental group decomposition for each component of the fuel

These files must have the same number of rows and the same order of components. Many examples can be found in the `fuelData <https://github.com/NREL/FuelLib/tree/main/fuelData>`_ directory.

Basic Usage
^^^^^^^^^^^

To demonstrate the usage of FuelLib, we will use the fuel "heptane-decane", which is a 
binary mixture of heptane and decane. The initial weight percentage composition is 73.75% 
heptane and 26.25% decane, and the group decomposition data is provided in the
`groupDecompositionData <https://github.com/NREL/FuelLib/tree/main/fuelData/groupDecompositionData>`_ directory.
The following tutorial is included in the `FuelLib/tutorials <https://github.com/NREL/FuelLib/tree/main/tutorials>`_
as ``basic.py``. To begin, we will import the necessary modules and create a ``groupContribution`` object for the two component fuel "heptane-decane": 

.. code-block:: python

    import os
    import sys

    # Add the FuelLib directory to the Python path
    FUELLIB_DIR = os.path.dirname(os.path.dirname(__file__))
    sys.path.append(FUELLIB_DIR)
    import paths
    import FuelLib as fl

    # Create a fuel object for the fuel "heptane-decane"
    fuel = fl.fuel("heptane-decane")

Upon initialization, the ``fuel`` object will read the initial weight 
percentage composition and group decomposition data from the specified files. The object stores
vectors of the calculated fundamental properties at standard conditions for each component of the fuel as described in :ref:`eq-GCM-properties`. 
For example, we can display the fuel name, the components in the fuel, the initial composition, and the critical temperature for each component: 

.. code-block:: python

    # Display fuel name, components, initial composition, and critical temperature
    print(f"Fuel name: {fuel.name}")
    print(f"Fuel components: {fuel.compounds}")
    print(f"Initial composition: {fuel.Y_0}")
    print(f"Critical temperature: {fuel.Tc} K")

.. code-block:: none

    >> Fuel name: heptane-decane
    >> Fuel components: ['NC7H16', 'NC10H22']
    >> Initial composition: [0.7375 0.2625]
    >> Critical temperature: [549.85598051 623.69051582] K

Next, we can calculate any of the component- or mixture-level properties using the 
``groupContribution`` object. For example, we can calculate the saturated vapor pressure
for each component and the mixture at a given temperature:

.. code-block:: python

    # Calculate the saturated vapor pressure at 320 K
    T = 320 # K
    p_sat_i = fuel.psat(T)
    p_sat_mix = fuel.mixture_vapor_pressure(T)
    print(f"Saturated vapor pressure at {T} K: {p_sat_i} Pa")
    print(f"Mixture saturated vapor pressure at {T} K: {p_sat_mix} Pa")

.. code-block:: none

    >> Saturated vapor pressure at 320 K: [13735.84605413   673.28876023] Pa
    >> Mixture saturated vapor pressure at 320 K: 11117.84926875165 Pa

The following links provide more information on the :ref:`eq-GCM-correlations` and
the :ref:`eq-mixture-properties` that can be calculated using the ``groupContribution`` object.

Exporting GCM Properties for Pele
---------------------------------

The development of FuelLib was motivated by the need for more accurate liquid fuel
property prediction in computational fluid dynamics (CFD) simulations. The fundamental GCM 
properties can be exported for use in the spray module of the `PelePhysics <https://github.com/AMReX-Combustion/PelePhysics>`_ library\ :footcite:p:`owen_pelemp_2024`
for combustion simulations in the `PeleLMeX <https://github.com/AMReX-Combustion/PeleLMeX>`_ 
flow solver\ :footcite:p:`henry_de_frahan_pele_2024` \ :footcite:p:`esclapez_pelelmex_2023`.

The export script, ``Export4Pele.py``, generates an input file named ``sprayPropsGCM_<fuel_name>.inp`` containing 
the necessary properties for each compound in the fuel. The properties are formatted for use in Pele and includes:

- Initial mass fraction
- Molecular weight
- Critical temperature
- Critical pressure
- Critical volume
- Boiling point
- Accentric factor
- Molar volume
- Specific heat
- Latent heat of vaporization

.. warning::
    The incorporation of the GCM in Pele is still under development and additional testing is required.

This example walks through the process and the available options for exporting GCM properties of a fuel named
"heptane-decane", which is a binary mixture of heptane and decane, using the ``Export4Pele.py`` script.

Default Options
^^^^^^^^^^^^^^^
.. note::
    The units for PeleLMeX are MKS while the units for PeleC are CGS. This is the same for 
    the spray inputs. Therefore, when running a spray simulation coupled with PeleC, the units for the 
    liquid fuel properties must be in CGS. The default units for the ``Export4Pele.py`` script is MKS, 
    but users can specify CGS by using the ``--units cgs`` option.
    
From the ``FuelLib`` directory, run the following command in the terminal, noting that ``--fuel_name`` is the only required input: ::
    
    cd FuelLib/source
    python Export4Pele.py --fuel_name heptane-decane


This generates the following input file, ``FuelLib/exportData/sprayPropsGCM_heptane-decane.inp``, for use in a PeleLMeX simulation: ::

    # -----------------------------------------------------------------------------
    # sprayPropsGCM_heptane-decane.inp
    # Generated on <YYY-MM-DD> <HH-MM-SS>
    # FuelLib remote URL: https://github.com/NREL/FuelLib.git
    # Git commit: <commit-hash>
    # Units: MKS
    # -----------------------------------------------------------------------------

    particles.spray_fuel_num = 2
    particles.fuel_species = NC7H16 NC10H22
    particles.Y_0 = 0.7375 0.2625
    particles.dep_fuel_names = NC7H16 NC10H22

    # Properties for NC7H16 in MKS
    particles.NC7H16_molar_weight = 0.100000 # kg/mol
    particles.NC7H16_crit_temp = 549.855981 # K
    particles.NC7H16_crit_press = 2821129.514417 # Pa
    particles.NC7H16_crit_vol = 0.000425 # m^3/mol
    particles.NC7H16_boil_temp = 379.073212 # K
    particles.NC7H16_acentric_factor = 0.336945 # -
    particles.NC7H16_molar_vol = 0.000146 # m^3/mol
    particles.NC7H16_cp_a = 1636.255000 # J/kg/K
    particles.NC7H16_cp_b = 3046.511000 # J/kg/K
    particles.NC7H16_cp_c = -983.629000 # J/kg/K
    particles.NC7H16_latent = 383110.000000 # J/kg

    # Properties for NC10H22 in MKS
    particles.NC10H22_molar_weight = 0.142000 # kg/mol
    particles.NC10H22_crit_temp = 623.690516 # K
    particles.NC10H22_crit_press = 2115522.932445 # Pa
    particles.NC10H22_crit_vol = 0.000592 # m^3/mol
    particles.NC10H22_boil_temp = 452.596977 # K
    particles.NC10H22_acentric_factor = 0.468050 # -
    particles.NC10H22_molar_vol = 0.000196 # m^3/mol
    particles.NC10H22_cp_a = 1630.488028 # J/kg/K
    particles.NC10H22_cp_b = 3098.105634 # J/kg/K
    particles.NC10H22_cp_c = -1024.456338 # J/kg/K
    particles.NC10H22_latent = 368035.211268 # J/kg

To include these parameters in your Pele simulation, copy the ``sprayPropsGCM_heptane-decane.inp`` 
file to the specific case directory and include the following line in your Pele input file: ::

    FILE = sprayPropsGCM.inp


Note: for liquid fuels from FuelLib with greater than 30 components, the script
will assume that all liquid fuel species deposit to the same gas-phase species, 
namely the name of the fuel. This is designed for conventional jet fuels such as POSF10325, where there are 
67 liquid fuel species corresponding to the GCxGC data, but only a single 
gas-phase mechanism species, "POSF10325". For example: ::

    cd FuelLib/source
    python Export4Pele.py --fuel_name posf10325

will result in the following: ::

    particles.spray_fuel_num = 67
    particles.fuel_species = Toluene C2-Benzene C3-Benzene ... C12-Tricycloparaffin
    particles.Y_0 = 0.001610 0.011172 0.0304982 ... 0.00110719
    particles.dep_fuel_names = POSF10325 POSF10325 ... POSF10325

    # Properties for Toluene in MKS
    ...

Additional Options
^^^^^^^^^^^^^^^^^^

There are four additional options that can be specified when running the export script:

- ``--units``: Specify the units for the properties. The default is "mks" but users can set the units to "cgs" for use in PeleC.
- ``--dep_fuel_names``: Specify which gas-phase species the liquid fuel deposits. The default is the same as the fuel name, but users can specify a single gas-phase species or a list of gas-phase species.
- ``--max_dep_fuels``: Specify the maximum number of dependent fuels. The default is 30 and is a bit arbitrary.
- ``--export_dir``: Specify the directory to export the file. The default is "FuelLib/exportData".

To specify all liquid fuel species deposity to a single gas-phase species, run the following command: ::

    cd FuelLib/source
    python Export4Pele.py --fuel_name heptane-decane --dep_fuel_names SINGLE_GAS

This will result in the following: ::

    particles.spray_fuel_num = 2
    particles.fuel_species = NC7H16 NC10H22
    particles.Y_0 = 0.7375 0.2625
    particles.dep_fuel_names = SINGLE_GAS SINGLE_GAS

    # Properties for NC7H16 in MKS
    ...

Alternatively, to specify a list of gas-phase species, run the following command: ::

    python Export4Pele.py --fuel_name heptane-decane --dep_fuel_names GAS_1 GAS_2

which produces: ::

    particles.spray_fuel_num = 2
    particles.fuel_species = NC7H16 NC10H22
    particles.Y_0 = 0.7375 0.2625
    particles.dep_fuel_names = GAS_1 GAS_2

    # Properties for NC7H16 in MKS
    ...

In the case that the liquid fuel has more than 30 components, the script will 
automatically set the deposition mapping to ``fuel.name`` for all components. 
If there are more than 30 components and the user wants each component to deposit 
to a gas-phase species of the same name, the user can increase ``--max_dep_fuels`` 
to a value greater than 30, however this would be required a massive mechanism for Pele and is not advised ::

    cd FuelLib/source
    python Export4Pele.py --fuel_name posf10325 --max_dep_fuels 67


Exporting GCM-Based Mixture Properties for Converge
---------------------------------------------------

The export script, ``Export4Converge.py``, generates a csv file named ``mixturePropsGCM_<fuel_name>.csv`` containing 
mixture property predictions for a given fuel over a specified temperature range. The properties include:

- Critical temperature
- Dynamic viscosity
- Surface tension
- Latent heat of vaporization
- Vapor pressure
- Density
- Specific heat
- Thermal conductivity

.. warning::
    Mixture properties for critical temperature, latent heat, and specific heat are provided by :ref:`conventional-mixing-rules` and need additional validation.

This example walks through the process and the available options for exporting GCM-based mixture properties for 
"posf10325", which is conventional Jet-A, using the ``Export4Converge.py`` script.

Default Options
^^^^^^^^^^^^^^^
    
From the ``FuelLib`` directory, run the following command in the terminal, noting that ``--fuel_name`` is the only required input: ::
    
    cd FuelLib/source
    python Export4Converge.py --fuel_name posf10325


This generates the file ``FuelLib/exportData/mixturePropsGCM_posf10325.csv`` with mixture 
property predictions from 0 K to 1000 K for use in a Converge simulation.

Additional Options
^^^^^^^^^^^^^^^^^^

There are four additional options that can be specified when running the export script:

- ``--units``: Specify the units for the mixture properties. The default is "mks" but users can set the units to "cgs".
- ``--temp_min``: Specify the minimum temperature. The default is 0 K.
- ``--temp_max``: Specify the maximum temperature. The default is 1000 K.
- ``--temp_step``: Specify the temperature step size. The default is :math:`\Delta T = 10` K.
- ``--export_dir``: Specify the directory to export the file. The default is "FuelLib/exportData".
  
.. note::
    The mixture property predictions may not be valid from the specified ``temp_min`` to ``temp_max``, 
    as the mixture properties are based on the GCM properties and correlations of the individual 
    components. Constant values are set for temperatures below the freezing point of the mixture or above 
    the minimum critical temperature of all compounds in the fuel. These temperature values will be noted in the 
    terminal output and should be considered when using the mixture properties in a simulation.


.. footbibliography::