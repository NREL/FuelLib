Tutorials
=========

This section provides tutorials on how to use the FuelLib library. 

Exporting GCM Properties for Pele
---------------------------------

The development of FuelLib was motivated by the need for more accurate liquid fuel
property prediction in computational fluid dynamics (CFD) simulations. The fundamental GCM 
properties can be exported for use in the spray module of the `PelePhysics <https://github.com/AMReX-Combustion/PelePhysics>`_ library\ :footcite:p:`owen_pelemp_2024`
for combustion simulations in the `PeleLMeX <https://github.com/AMReX-Combustion/PeleLMeX>`_ 
flow solver\ :footcite:p:`henry_de_frahan_pele_2024` \ :footcite:p:`esclapez_pelelmex_2023`.

The export script, ``Export4Pele.py``, generates an input file named ``sprayPropsGCM.inp`` containing 
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
From the ``FuelLib`` directory, run the following command in the terminal, noting that ``--fuel_name`` is the only required input: ::
    
    python Export4Pele.py --fuel_name heptane-decane


This generates the following input file ``FuelLib/sprayPropsGCM/sprayPropsGCM.inp``: ::

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
    particles.NC7H16_cp = 1636.255 3046.5109999999995 -983.6289999999999 # J/kg/K
    particles.NC7H16_latent = 383110.000000 # J/kg

    # Properties for NC10H22 in MKS
    particles.NC10H22_molar_weight = 0.142000 # kg/mol
    particles.NC10H22_crit_temp = 623.690516 # K
    particles.NC10H22_crit_press = 2115522.932445 # Pa
    particles.NC10H22_crit_vol = 0.000592 # m^3/mol
    particles.NC10H22_boil_temp = 452.596977 # K
    particles.NC10H22_acentric_factor = 0.468050 # -
    particles.NC10H22_molar_vol = 0.000196 # m^3/mol
    particles.NC10H22_cp = 1630.488028169014 3098.1056338028166 -1024.456338028169 # J/kg/K
    particles.NC10H22_latent = 368035.211268 # J/kg

To include these parameters in your Pele simulation, copy the ``sprayPropsGCM.inp`` 
file to the specific case directory and include the following line in your Pele input file: ::

    FILE = sprayPropsGCM.inp


Note: for liquid fuels from FuelLib with greater than 30 components, the script
will assume that all liquid fuel species deposit to the same gas-phase species, 
namely the name of the fuel. This is designed for conventional jet fuels such as POSF10325, where there are 
67 liquid fuel species correpsonding to the GCxGC data, but only a single 
gas-phase mechanism species, "POSF10325". For example: ::

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
- ``--export_dir``: Specify the directory to export the file. The default is "FuelLib/sprayPropsGCM".

To specify all liquid fuel species deposity to a single gas-phase species, run the following command: ::

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

    python Export4Pele.py --fuel_name posf10325 --max_dep_fuels 67


.. footbibliography::