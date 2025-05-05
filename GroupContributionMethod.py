import pandas as pd
import numpy as np
import os


class groupContribution:
    """
    Class for handling group contribution calculations of thermodynamic properties
    and mixture properties.
    """

    # Paths to input directories
    FuelLib_PATH = os.path.dirname(os.path.abspath(__file__))
    gcmTableDir = os.path.join(FuelLib_PATH, "gcmTableData")
    fuelDataDir = os.path.join(FuelLib_PATH, "fuelData")
    groupDecompDir = os.path.join(fuelDataDir, "groupDecompositionData")
    gcxgcDir = os.path.join(fuelDataDir, "gcData")

    # Path to GCM table
    gcmTableFile = os.path.join(gcmTableDir, "gcmTable.csv")

    # Class-level variables to hold mixture-specific data and GCM table properties
    name = ""
    num_compounds = None

    # Initial composition and functional group data for mixture
    Y_0 = None  # Initial mass fraction for mixture (num_compounds,)
    Nij = None  # Compound vs. group matrix (num_compounds, num_groups)

    # Number of first and second order groups from Constantinou and Gani
    N_g1 = 78
    N_g2 = 43

    # critical properties at standard temp (num_compounds,)
    MW = None
    Tc = None
    Pc = None
    Vc = None
    Tb = None
    Tm = None
    Hf = None
    Gf = None
    Hv_stp = None
    omega = None
    Vm_stp = None
    Cp_stp = None
    Cp_B = None
    Cp_C = None
    Lv_stp = None

    # Classifying hydrocarbon family
    fam = None

    # Boltzmann's constant J/K
    k_B = 1.380649e-23

    def __init__(self, name, decompName=None, W=1):
        """
        Initializes mixture-specific data and GCM properties for the specified
        mixture. Reads GCM table and mixture data from files.

        Parameters:
        name (str): Name of the mixture as it appears in its gcData file.
        decompName (str, optional): Name of the groupDecomposition file if different from name.
        W (int, optional): Determines if first-order only approximation (i.e. W = 0)
        """

        self.name = name
        if decompName is None:
            decompName = name
        groupDecompFile = os.path.join(self.groupDecompDir, f"{decompName}.csv")
        gcxgcFile = os.path.join(self.gcxgcDir, f"{name}_init.csv")

        # Read functional group data for mixture (num_compounds,num_groups)
        df_Nij = pd.read_csv(groupDecompFile)
        self.Nij = df_Nij.iloc[:, 1:].to_numpy()
        self.num_compounds = self.Nij.shape[0]
        self.num_groups = self.Nij.shape[1]

        # Classify hydrocarbon by family (used in thermal conductivity)
        # 0: saturated hydrocarbons
        # 1: aromatics
        # 2: cycloparaffins
        # 3: olefins
        self.fam = np.zeros(self.num_compounds, dtype=int)
        aromatics = 10  # starting index for aromatic groups
        num_aromatics = 5
        cyclos = 84  # starting index for membered ring groups
        num_cyclos = 5
        olefins = 4  # starting index for double bound groups
        num_olefins = 6
        for i in range(self.num_compounds):
            # Check if aromatic: does it contain AC's?
            if sum(self.Nij[i, aromatics : aromatics + num_aromatics]) > 0:
                self.fam[i] = 1
            # Check if cycloparaffin: does it contain rings?
            elif sum(self.Nij[i, cyclos : cyclos + num_cyclos]) > 0:
                self.fam[i] = 2
            # Check if olefin: does it contain double bonds?
            elif sum(self.Nij[i, olefins : olefins + num_olefins]) > 0:
                self.fam[i] = 3

        # Read initial liquid composition of mixture and normalize to get mass frac
        df_gcxgc = pd.read_csv(gcxgcFile, usecols=[1])
        self.Y_0 = df_gcxgc.to_numpy().flatten().astype(float)
        self.Y_0 /= np.sum(self.Y_0)

        # Make sure mixture data is consistent:
        if self.num_groups < self.N_g1:
            raise ValueError(
                f"Insufficient mixture description:\n"
                f"The number of columns in {groupDecompFile} is less than "
                f"the required number of first-order groups (N_g1 = {self.N_g1})."
            )
        if self.Y_0.shape[0] != self.num_compounds:
            raise ValueError(
                f"Insufficient mixture description:\n"
                f"The number of compounds in {groupDecompFile} does not "
                f"equal the number of compounds in {gcxgcFile}."
            )

        # Read and store GCM table properties
        df_gcm_properties = pd.read_csv(self.gcmTableFile)
        gcm_properties = df_gcm_properties.loc[
            :, ~df_gcm_properties.columns.isin(["Property", "Units"])
        ].to_numpy()

        # Determine if approximations include second-order contributions
        if W == 0 or self.num_groups < self.N_g1 + self.N_g2:
            # Only retain first-order GCM properties
            gcm_properties = gcm_properties[:, 0 : self.N_g1]
            self.Nij = self.Nij[:, 0 : self.N_g1]

        # Table data for functional groups (num_compounds,)
        Tck = gcm_properties[0]  # critical temperature (1)
        Pck = gcm_properties[1]  # critical pressure (bar)
        Vck = gcm_properties[2]  # critical volume (m^3/kmol)
        Tbk = gcm_properties[3]  # boiling temperature (1)
        Tmk = gcm_properties[4]  # melting point temperature (1)
        hfk = gcm_properties[5]  # enthalpy of formation, (kJ/mol)
        gfk = gcm_properties[6]  # Gibbs energy (kJ/mol)
        hvk = gcm_properties[7]  # latent heat of vaporization (kJ/mol)
        wk = gcm_properties[8]  # accentric factor (1)
        Vmk = gcm_properties[9]  # liquid molar volume fraction (m^3/kmol)
        cpak = gcm_properties[10]  # specific heat values (J/mol/K)
        cpbk = gcm_properties[11]  # specific heat values (J/mol/K)
        cpck = gcm_properties[12]  # specific heat values (J/mol/K)
        mwk = gcm_properties[13]  # molecular weights (g/mol)

        # --- Compute critical properties at standard temp (num_compounds,)
        # Molecular weights
        self.MW = np.matmul(self.Nij, mwk)  # g/mol
        self.MW *= 1e-3  # Convert to kg/mol

        # T_c (critical temperature)
        self.Tc = 181.128 * np.log(np.matmul(self.Nij, Tck))  # K

        # p_c (critical pressure)
        self.Pc = 1.3705 + (np.matmul(self.Nij, Pck) + 0.10022) ** (-2)  # bar
        self.Pc *= 1e5  # Convert to Pa from bar

        # V_c (critical volume)
        self.Vc = -0.00435 + (np.matmul(self.Nij, Vck))  # m^3/kmol
        self.Vc *= 1e-3  # Convert to m^3/mol

        # T_b (boiling temperature)
        self.Tb = 204.359 * np.log(np.matmul(self.Nij, Tbk))  # K

        # T_m (melting temperature)
        self.Tm = 102.425 * np.log(np.matmul(self.Nij, Tmk))  # K

        # H_f (enthalpy of formation)
        self.Hf = 10.835 + np.matmul(self.Nij, hfk)  # kJ/mol
        self.Hf *= 1e3  # Convert to J/mol

        # G_f (Gibbs free energy)
        self.Gf = -14.828 + np.matmul(self.Nij, gfk)  # kJ/mol
        self.Gf *= 1e3  # Convert to J/mol

        # H_v,stp (enthalpy of vaporization at 298 K)
        self.Hv_stp = 6.829 + (np.matmul(self.Nij, hvk))  # kJ/mol
        self.Hv_stp *= 1e3  # Convert to J/mol

        # omega (accentric factor)
        self.omega = 0.4085 * np.log(np.matmul(self.Nij, wk) + 1.1507) ** (1.0 / 0.5050)

        # V_m (molar liquid volume at 298 K)
        self.Vm_stp = 0.01211 + np.matmul(self.Nij, Vmk)  # m^3/kmol
        self.Vm_stp *= 1e-3  # Convert to m^3/mol

        # C_p,stp (specific heat at 298 K)
        self.Cp_stp = np.matmul(self.Nij, cpak) - 19.7779  # J/mol/K

        # Temperature corrections for C_p
        self.Cp_B = np.matmul(self.Nij, cpbk)
        self.Cp_C = np.matmul(self.Nij, cpck)

        # L_v,stp (latent heat of vaporization at 298 K)
        self.Lv_stp = self.Hv_stp / self.MW  # J/kg

        # Lennard-Jones parameters for diffusion calculations (Tee et al. 1966)
        self.epsilonByKB = (0.7915 + 0.1693 * self.omega) * self.Tc  # K
        Pc_atm = self.Pc / 101325  # atm
        self.sigma = (2.3551 - 0.0874 * self.omega) * (self.Tc / Pc_atm) ** (
            1.0 / 3
        )  # Angstroms
        self.sigma *= 1e-10  # Convert from Angstroms to m

    # -------------------------------------------------------------------------
    # Member functions
    # -------------------------------------------------------------------------
    def mass_frac(self, mass):
        """
        Calculate the mass fractions from the mass of each component.

        Parameters:
        mass (np.ndarray): Mass of each compound (shape: num_compounds,).

        Returns:
        Yi np.ndarray: Mass fractions of the compounds (shape: num_compounds,).
        """
        # Normalize to get group mole fractions
        total_mass = np.sum(mass)
        if total_mass != 0:
            Yi = mass / total_mass
        else:
            Yi = np.zeros_like(self.MW)

        return Yi

    def mole_frac(self, mass):
        """
        Calculate the mole fractions from the mass of each component.

        Parameters:
        mass (np.ndarray): Mass of each compound in kg (shape: num_compounds,).

        Returns:
        Xi np.ndarray: Molar fractions of the compounds (shape: num_compounds,).
        """
        # Calculate the number of moles for each compound
        num_mole = mass / self.MW

        # Normalize to get group mole fractions
        total_moles = np.sum(num_mole)
        if total_moles != 0:
            Xi = num_mole / total_moles
        else:
            Xi = np.zeros_like(self.MW)

        return Xi

    def density(self, T):
        """
        Calculate the density of each component at temperature T.

        Parameters:
        T (float): Temperature of the mixture in Kelvin

        Return:
        rho (np.array): Density of each compound in kg/m^3 (shape: num_compounds,)
        """
        MW = self.MW  # kg/mol
        Vm = self.molar_liquid_vol(T)  # m^3/mol
        rho = MW / Vm  # Density of each component of the mixture kg/m^3

        return rho

    def viscosity_kinematic(self, T):
        """
        Calculate the viscosity of individual components in the mixture at a given
        temperature using Dutt's equation (4.23) in Viscosity of Liquids.

        Parameters:
        T (float): Temperature in Kelvin.

        Returns:
        np.array: Viscosity of each component in m^2/s.
        """

        # Convert temperature to Celsius
        T_cels = K2C(T)
        Tb_cels = K2C(self.Tb)

        # RHS of Dutt's equation (4.23) in Viscosity of Liquids
        rhs = -3.0171 + (442.78 + 1.6452 * Tb_cels) / (T_cels + 239 - 0.19 * Tb_cels)
        nu_i = np.exp(rhs)  # Viscosity in mm^2/s

        # Convert to SI (m^2/s)
        nu_i *= 1e-6

        return nu_i

    def viscosity_dynamic(self, T):
        """
        Calculates liquid dynamic viscosity based on droplet temperature and
        density using Dutt's Equation (4.23) in "Viscosity of Liquids".

        Parameters:
        T (float): Temperature in Kelvin.

        Returns:
        mu (np.ndarray): Dynamic viscosity in Pa*s (shape: num_compounds,).
        """
        nu_i = self.viscosity_kinematic(T)  # m^2/s
        rho_i = self.density(T)  # kg/m^3
        mu_i = nu_i * rho_i  # Pa*s
        return mu_i

    def Cp(self, T):
        """
        Computes specific heat capacity at a given temperature.

        Parameters:
        T (float): Temperature in Kelvin.

        Returns:
        np.ndarray: Specific heat capacity in J/mol/K (shape: num_compounds,).
        """
        theta = (T - 298) / 700
        cp = self.Cp_stp + self.Cp_B * theta + self.Cp_C * theta**2
        return cp

    def Cl(self, T):
        """
        Computes liquid specific heat capacity at a given temperature.

        Parameters:
        T (float): Temperature in Kelvin.

        Returns:
        np.ndarray: Specific heat capacity in J/kg/K (shape: num_compounds,).
        """
        cp = self.Cp(T)
        return cp / self.MW

    def psat(self, T, correlation="Lee-Kesler"):
        """
        Computes the saturated vapor pressure.

        Parameters:
        T (float): Temperature in Kelvin.
        correlation (str, optional): "Ambrose-Walton" or "Lee-Kesler".

        Returns:
        np.ndarray: Saturated vapor pressure in Pa.
        """
        Tr = T / self.Tc

        if correlation.casefold() == "Ambrose-Walton".casefold():
            # May cause trouble at high temperatures
            tau = 1 - Tr
            f0 = (
                -5.97616 * tau
                + 1.29874 * tau**1.5
                - 0.60394 * tau**2.5
                - 1.06841 * tau**5.0
            )
            f0 /= Tr
            f1 = (
                -5.03365 * tau
                + 1.11505 * tau**1.5
                - 5.41217 * tau**2.5
                - 7.46628 * tau**5.0
            )
            f1 /= Tr
            f2 = (
                -0.64771 * tau
                + 2.41539 * tau**1.5
                - 4.26979 * tau**2.5
                - 3.25259 * tau**5.0
            )
            f2 /= Tr
            rhs = np.exp(f0 + self.omega * f1 + self.omega**2 * f2)

        else:  # Default correlation is Lee-Kesler
            f0 = 5.92714 - (6.09648 / Tr) - 1.28862 * np.log(Tr) + 0.169347 * (Tr**6)
            f1 = 15.2518 - (15.6875 / Tr) - 13.4721 * np.log(Tr) + 0.43577 * (Tr**6)
            rhs = np.exp(f0 + self.omega * f1)

        psat = self.Pc * rhs
        return psat

    def molar_liquid_vol(self, T):
        """
        Computes molar liquid volumes with temperature correction.

        Parameters:
        T (float): Temperature in Kelvin.

        Returns:
        np.ndarray: Molar liquid volume in m^3/mol (num_compounds,)
        """

        Tstp = 298.0
        phi = np.zeros_like(self.Tc)
        for i in range(len(self.Tc)):
            if T > self.Tc[i]:
                phi[i] = -((1 - (Tstp / self.Tc[i])) ** (2.0 / 7.0))
            else:
                phi[i] = ((1 - (T / self.Tc[i])) ** (2.0 / 7.0)) - (
                    (1 - (Tstp / self.Tc[i])) ** (2.0 / 7.0)
                )
        z_vec = 0.29056 - 0.08775 * self.omega
        Vmi = self.Vm_stp * np.power(z_vec, phi)
        return Vmi

    def latent_heat_vaporization(self, T):
        """
        Calculates the letent heat of vaporization, adjusted for temperature.

        Parameters:
        T (float): Temperature in Kelvin.

        Returns:
        np.ndarray: Adjusted latent heat of vaporization in J/kg (num_compounds,)
        """
        Lvi = np.zeros_like(self.Tc)

        Tr = T / self.Tc
        Trb = self.Tb / self.Tc

        for i in range(0, self.Tc.shape[0]):
            if T > self.Tc[i]:
                Lvi[i] = 0.0
            else:
                Lvi[i] = self.Lv_stp[i] * (((1.0 - Tr[i]) / (1.0 - Trb[i])) ** 0.38)

        return Lvi

    def diffusion_coeff(
        self,
        p,
        T,
        sigma_gas=3.62e-10,
        epsilonByKB_gas=97.0,
        MW_gas=28.97e-3,
        correlation="Tee",
    ):
        """
        Computes the diffusion coefficient using Lennard-Jones parameters according
        to Wilke and Lee.  See equation (11-4.1) in Poling. Assumes ambient gas is air
        with values from Poling (see end of Section 11-4.1).

        Parameters:
        p (float): Pressure in Pa.
        T (float): Temperature in Kelvin.
        sigma_gas (float, optional): collision diameter (m), default is air
        epsilonByKB_gas (float, optional): (well depth) / (Boltzmann_const) (K), default is air
        MW_gas (float, optional): Mean molecular weight of gas (kg/mol), default is air
        correlation (str, optional): Method for calculating sigma and epsilonByKB

        Returns:
        D_AB (np.ndarray): Diffusion coefficient (num_compounds,).
        """

        # Method of Tee for calculating liquid sigma and epsilon
        if correlation.casefold() == "Tee".casefold():
            sigma_i = self.sigma * 1e10  # convert from m to Angstroms
            epsilonByKB_i = self.epsilonByKB  # K
        else:
            # Method of Wilke & Lee calculating liquid sigma and epsilon
            Vmb_i = np.zeros_like(self.Tb)
            for n in range(self.num_compounds):
                Vmb_i[n] = self.molar_liquid_vol(self.Tb[n])[n] * 1e6  # cm^3/mol
            sigma_i = 1.18 * Vmb_i ** (1 / 3)  # Angstroms, Poling (11-4.2)
            epsilonByKB_i = 1.15 * self.Tb  # K , Poling (11-4.3)

        # Compute binary sigma and epsilon
        sigma_gas *= 1e10  # convert from m to Angstroms
        sigmaAB_i = (sigma_gas + sigma_i) / 2  # Angstroms, Poling (11-3.5)
        epsilonAB_byKB_i = (
            epsilonByKB_gas * epsilonByKB_i
        ) ** 0.5  # K, Poling (11-3.4)

        # Dimensionless collision integral for diffusion: Poling (11-3.6)
        Tstar_i = T / epsilonAB_byKB_i  # [1]
        A = 1.06036
        B = 0.15610
        C = 0.193
        D = 0.47635
        E = 1.03587
        F = 1.52996
        G = 1.76474
        H = 3.89411
        omegaD_i = (
            A / (Tstar_i**B)
            + C / np.exp(D * Tstar_i)
            + E / np.exp(F * Tstar_i)
            + G / np.exp(H * Tstar_i)
        )

        # Convert molecular weights from kg/mol to g/mol then calculate M_AB
        MW_gas *= 1e3
        MW_i = self.MW * 1e3
        M_AB_i = 2 * (MW_i * MW_gas) / (MW_i + MW_gas)  # g/mol, see Poling (11-3.1)

        # Convert pressure from Pa to bar
        p *= 1e-5  # bar

        # Binary diffusion coefficients, Poling (11-4.1)
        D_AB_i = (
            1e-3
            * (3.03 - 0.98 / (M_AB_i**0.5))
            * (T**1.5)
            / (p * M_AB_i**0.5 * sigmaAB_i**2 * omegaD_i)
        )  # cm^2/s
        D_AB_i *= 1e-4  # Convert to m^2/s

        return D_AB_i

    def surface_tension(self, T, correlation="Brock-Bird"):
        """
        Calculate the surface tension at a given temperature. Both
        correlations can be found in Poling's book. Brock-Bird is eqt. 12-3.5 and
        Pitzer is eqt. 12-3.7.

        Parameters:
        T (float): Temperature in Kelvin.
        correlation (str, optional): "Brock-Bird" or "Pitzer".

        Returns:
        np.ndarray: Surface tension in N/m (num_compounds,)
        """
        Tr = T / self.Tc
        Pc = self.Pc * 1e-5  # convert from Pa to bar
        Tc = self.Tc  # K

        if correlation.casefold() == "Brock-Bird".casefold():
            Tbr = self.Tb / Tc
            Q = 0.1196 * (1.0 + (Tbr * np.log(Pc / 1.01325)) / (1.0 - Tbr)) - 0.279
        else:
            w = self.omega
            Q = (
                (1.86 + 1.18 * w)
                / 19.05
                * (((3.75 + 0.91 * w) / (0.291 - 0.08 * w)) ** (2.0 / 3.0))
            )

        st = Pc ** (2.0 / 3.0) * Tc ** (1.0 / 3.0) * Q * (1 - Tr) ** (11.0 / 9.0)

        st *= 1e-3  # Convert from dyn/cm to N/m

        return st

    def thermal_conductivity(self, T):
        """
        Calculate the thermal conductivity at a given temperature using the
        method of Latini et al. equation (10-9.1) in Poling.

        Parameters:
        T (float): Temperature in Kelvin.

        Returns:
        np.ndarray: Thermal conductivity in W/m/K (num_compounds,)
        """

        Astar = 0.00350 + np.zeros(self.num_compounds)
        alpha = 1.2
        beta = 0.5 + np.zeros(self.num_compounds)
        gamma = 0.167
        MW_beta = self.MW * 1e3  # convert from kg/mol to g/mol
        Tr = T / self.Tc

        for i in range(self.num_compounds):
            if self.fam[i] == 1:
                # Aromatics
                Astar[i] = 0.0346
                beta[i] = 1.0
            elif self.fam[i] == 2:
                # Cycloparaffins
                Astar[i] = 0.0310
                beta[i] = 1.0
            elif self.fam[i] == 3:
                # Olefins
                Astar[i] = 0.0361
                beta[i] = 1.0
            MW_beta[i] = MW_beta[i] ** beta[i]

        A = Astar * self.Tb**alpha / (MW_beta * self.Tc**gamma)
        tc = A * (1 - Tr) ** (0.38) / (Tr ** (1 / 6))

        return tc

    # --- Mixture functions ---
    def mixture_density(self, Yi, T):
        """
        Calculate the mixture density at a given temperature.

        Parameters:
        Yi (np.ndarray): Mass fractions of each compound (shape: num_compounds,).
        T (float): Temperature in Kelvin.

        Returns:
        float: Mixture density in kg/m^3.
        """
        MW = self.MW  # Molecular weights of each component (kg/mol)
        Vmi = self.molar_liquid_vol(T)  # Molar volume of each component (m^3/mol)

        # Calculate density (kg/m^3)
        rho = Yi @ (MW / Vmi)

        return rho

    def mixture_kinematic_viscosity(self, mass, T, correlation="Kendall-Monroe"):
        """
        Calculate the kinematic viscosity of the mixture at a given temperature.

        Parameters:
        mass (np.ndarray): Mass of each compound of mixture (shape: num_compounds,).
        T (float): Temperature in Kelvin.
        correlation (str, optional): Mixing model "Kendall-Monroe" or "Arrhenius".

        Returns:
        float: Mixture viscosity in mm^2/s.
        """
        nu_i = self.viscosity_kinematic(T)  # Viscosities of individual components

        # Calculate group mole fractions for each species
        Xi = self.mole_frac(mass)

        if correlation.casefold() == "Arrhenius".casefold():
            # Arrhenius mixing correlation
            nu = np.exp(np.sum(Xi * np.log(nu_i)))
        else:
            # Default: Kendall-Monroe mixing correlation
            nu = np.sum(Xi * (nu_i ** (1.0 / 3.0))) ** (3.0)

        return nu

    def mixture_dynamic_viscosity(self, mass, T, correlation="Kendall-Monroe"):
        """
        Calculate the dynamic viscosity of the mixture at a given temperature.

        Parameters:
        mass (np.ndarray): Mass of each compound of mixture (shape: num_compounds,).
        T (float): Temperature in Kelvin.
        correlation (str, optional): Mixing model "Kendall-Monroe", "Arrhenius".

        Returns:
        float: Mixture viscosity in Pa*s.
        """
        nu = self.mixture_kinematic_viscosity(mass, T, correlation)
        Yi = self.mass_frac(mass)
        rho = self.mixture_density(Yi, T)

        return rho * nu

    def mixture_vapor_pressure(self, mass, T, correlation="Lee-Kesler"):
        """
        Calculate the vapor pressure the mixture.

        Parameters:
        mass (np.ndarray): Mass of each compound of mixture (shape: num_compounds,).
        T (float): Temperature in Kelvin.
        correlation (str, optional): "Ambrose-Walton" or "Lee-Kesler".

        Returns:
        float: vapor pressure in Pa.
        """

        # Group mole fraction for each compound
        Xi = self.mole_frac(mass)

        # Saturated vapor pressure for each compound (Pa)
        p_sati = self.psat(T, correlation)

        # Mixture vapor pressure via Raoult's law
        p_v = p_sati @ Xi

        return p_v

    def mixture_surface_tension(self, mass, T, correlation="Brock-Bird"):
        """
        Calculate the surface tension of the mixture using an arithmetic psuedo-
        property as recommended by Hugill and van Welsenes (1986).

        Parameters:
        mass (np.ndarray): Mass of each compound of mixture (shape: num_compounds,).
        T (float): Temperature in Kelvin.
        correlation (str, optional): "Pitzer" or "Brock-Bird".

        Returns:
        float: surface tension in N/m.
        """

        # Group mole fraction for each compound
        Xi = self.mole_frac(mass)

        # Surface tension for each compound (N/m)
        sti = self.surface_tension(T, correlation)

        # Mixture surface tension via arithmetic mean, Poling (12-5.2)
        st = mixingRule(sti, Xi, "arithmetic")

        return st

    def mixture_thermal_conductivity(self, mass, T):
        """
        Calculate the thermal conductivity of the mixture.

        Parameters:
        T (float): Temperature in Kelvin.

        Returns:
        float: thermal conductivity in W/m/K.
        """
        Yi = self.mass_frac(mass)
        tc = self.thermal_conductivity(T)
        return np.sum(Yi * tc ** (-2)) ** (-0.5)


# -----------------------------------------------------------------------------
# Utility functions
# -----------------------------------------------------------------------------
def C2K(T):
    """
    Convert temperature from Celsius to Kelvin.

    Parameters:
    T (float or np.array): Temperature in Celsius.

    Returns:
    float or np.array: Temperature in Kelvin.
    """
    return T + 273.15


def K2C(T):
    """
    Convert temperature from Kelvin to Celsius.

    Parameters:
    T (float or np.array): Temperature in Kelvin.

    Returns:
    float or np.array: Temperature in Celsius.
    """
    return T - 273.15


def mixingRule(var_n, X, pseudo_prop="arithmetic"):
    """
    Mixing rules for computing mixture properties from individual compound properties.

    Parameters:
    var_n (np.array): Individual compound properties (shape: num_compounds,).
    X (np.array): Mole fraction (shape: num_compounds,).
    pseudo_prop (str, optional): Type of mean (e.g. arithmetic of geometric).

    Returns:
    float: Single mixture property with units of var_n.
    """
    num_comps = len(var_n)
    var_mix = 0.0
    for i in range(num_comps):
        for j in range(num_comps):
            if pseudo_prop.casefold() == "geometric":
                # Use geometric mean definition for the pseudo property
                var_ij = (var_n[i] * var_n[j]) ** (0.5)
            else:
                # Use arithmetic definition for the pseudo property
                var_ij = (var_n[i] + var_n[j]) / 2
            var_mix += X[i] * X[j] * var_ij
    return var_mix
