import numpy as np
from scipy.integrate import odeint
from config import Config


class AptamerBindingModel:
    """
    Encapsulates the math/physics/chemistry for an aptamer-based sensor:
    - Binding kinetics
    - Hill-Langmuir model
    - SWV frequency response
    - Butler-Volmer equation
    - Marcus theory (in eV)
    """
    
    F = Config.FARADAY_CONST                                            # Faraday's constant (C/mol)
    R = Config.UNI_GAS_CONST                                            # Universal gas constant (J/(mol·K))
    kB = Config.BOLTZMANN_JK_CONST                                      # Boltzmann constant (J/K)

    def __init__(self):
        # --------------------------------------------------------------------------------------------------------------------------------
        # 1. Binding Kinetics Parameters
        # --------------------------------------------------------------------------------------------------------------------------------
        self.A_total = Config.BIND_KIN_APT_CONC_TOTAL                   # Total aptamer concentration (arbitrary units)
        self.T_conc  = Config.BIND_KIN_TRGT_CONC                        # Target concentration
        self.k_on    = Config.BIND_KIN_K_ON                             # Binding rate constant (1/s)
        self.k_off   = Config.BIND_KIN_K_OFF                            # Unbinding rate constant (1/s)
        self.I_min   = Config.BIND_KIN_I_MIN                            # Minimum current
        self.I_max   = Config.BIND_KIN_I_MAX                            # Maximum current
        self.AT0     = Config.BIND_KIN_APT_INIT_CONC                    # Initial bound aptamer concentration
        self.n       = Config.BIND_KIN_HILL_COEFF                       # Hill coefficient

        # --------------------------------------------------------------------------------------------------------------------------------
        # 2. SWV Simulation Parameters
        # --------------------------------------------------------------------------------------------------------------------------------
        self.k_et_unbound = Config.SWV_SIM_K_ET_UNBOUND
        self.k_et_bound   = Config.SWV_SIM_K_ET_BOUND
        self.f_decay      = Config.SWV_SIM_FREQ_DECAY
        self.I0           = Config.SWV_SIM_I0

        # --------------------------------------------------------------------------------------------------------------------------------
        # 3. Butler–Volmer Parameters
        # --------------------------------------------------------------------------------------------------------------------------------
        self.j0    = Config.BV_CURR_DENS_J0                             # Exchange current density (A/cm²)
        self.alpha = Config.BV_ALPHA                                    # Charge transfer coefficient
        self.T     = Config.BV_TEMP                                     # Temperature (K)
        self.eta0  = Config.BV_ETA_0                                    # Baseline overpotential (V) when no target is bound.

        # --------------------------------------------------------------------------------------------------------------------------------
        # 4. Marcus Theory Parameters (in eV)
        # --------------------------------------------------------------------------------------------------------------------------------
        self.A_marcus      = Config.MARCUS_PREEXP_A                     # Pre-exponential factor (s⁻¹)
        self.lambda_marcus = Config.MARCUS_REORG_E_LAMBDA               # Reorganisation energy (eV)
        self.deltaG0       = Config.MARCUS_DELTAG                       # Standard free energy change (eV)
        self.deltaG_shift  = Config.MARCUS_DELTAG_SHIFT                 # Shift in free energy upon binding (eV)

        # Experimental data
        self.freq_low = Config.EXPRM_DATA_FREQ_LOW                      # Frequency below the critical frequency (e.g. 16 Hz)
        self.freq_high = Config.EXPRM_DATA_FREQ_HIGH                    # Frequency above the critical frequency (e.g. 240 Hz)
        
        self.charge_low_unbound  = Config.EXPRM_DATA_CHRG_LOW_UNBOUND   # Calculated charge (from integration) at these frequencies (in, say, µC)
        self.charge_low_bound    = Config.EXPRM_DATA_CHRG_LOW_BOUND
        self.charge_high_unbound = Config.EXPRM_DATA_CHRG_HIGH_UNBOUND
        self.charge_high_bound   = Config.EXPRM_DATA_CHRG_HIGH_BOUND
        
        self.electrode_area = Config.EXPRM_DATA_ELEC_AREA               # Electrode geometry (nominal area in mm²)
        self.mb_electron = Config.EXPRM_DATA_MB_N_ET                    # Number of electrons transferred during MB redox (default 2)


    ###########################################################################
    # 1. Aptamer Binding Kinetics
    ###########################################################################
    def _binding_ode(self, AT, t):
        A_free = self.A_total - AT
        return self.k_on * A_free * self.T_conc - self.k_off * AT

    def compute_binding_curve(self, 
                                t_end=Config.BIND_KIN_T_END,
                                num_points=Config.BIND_KIN_NUM_PTS
                                ):
        t = np.linspace(0, t_end, num_points)
        AT = odeint(self._binding_ode, self.AT0, t).flatten()
        fraction_bound = AT / (self.A_total + Config.EPSILON_CONST)
        I = self.I_min + (self.I_max - self.I_min) * fraction_bound
        return t, fraction_bound, I

    ###########################################################################
    # 2. Hill–Langmuir Isotherm
    ###########################################################################
    def hill_langmuir(self, T, Kd, n):
        return np.power(T, n) / (np.power(Kd, n) + np.power(T, n) + Config.EPSILON_CONST)

    def compute_hill_langmuir_curve(self,
                                    t_min=Config.HL_ISO_T_MIN,
                                    t_max=Config.HL_ISO_T_MAX,
                                    num_points=Config.HL_ISO_NUM_PTS
                                    ):
        Kd = self.k_off / (self.k_on + Config.EPSILON_CONST)
        T_conc_range = np.logspace(np.log10(t_min), np.log10(t_max), num_points)
        fraction_bound = self.hill_langmuir(T_conc_range, Kd, self.n)
        I_hl = self.I_min + (self.I_max - self.I_min) * fraction_bound
        return T_conc_range, fraction_bound, I_hl

    ###########################################################################
    # 3. Square Wave Voltammetry (SWV) - inverse frequency transform
    ###########################################################################

    def compute_swv_curve(self, num_points=Config.SWV_PARAM_NUM_PTS):
        # Use the experimental parameters for frequency bounds
        frequencies = np.linspace(self.freq_low, self.freq_high, num_points)
        
        # Calculate net SWV current (for unbound vs. bound)
        denom_unbound = (2 * frequencies + Config.EPSILON_CONST)
        denom_bound   = (2 * frequencies + Config.EPSILON_CONST)
        I_net_unbound = self.I0 * np.exp(-self.k_et_unbound / denom_unbound) * np.exp(-frequencies / (self.f_decay + Config.EPSILON_CONST))
        I_net_bound   = self.I0 * np.exp(-self.k_et_bound   / denom_bound)   * np.exp(-frequencies / (self.f_decay + Config.EPSILON_CONST))
        
        # Transform: (peak current) / frequency
        I_u_F = I_net_unbound / frequencies
        I_b_F = I_net_bound   / frequencies
        
        # x-axis: 1/frequency
        inv_freq = 1 / frequencies
        
        # Normalize for plotting
        I_u_norm = I_u_F / (np.max(I_u_F) + Config.EPSILON_CONST)
        I_b_norm = I_b_F / (np.max(I_b_F) + Config.EPSILON_CONST)
        
        # Find crossover point
        diff_array = np.abs(I_u_norm - I_b_norm)
        crossover_index = np.argmin(diff_array)
        crossover_freq = frequencies[crossover_index]
        crossover_inv_freq = inv_freq[crossover_index]
        
        return inv_freq, I_u_F, I_b_F, I_u_norm, I_b_norm, crossover_freq, crossover_inv_freq


    ###########################################################################
    # 4. Butler–Volmer Equation
    ###########################################################################
    def butler_volmer(self, eta):
        term_forward = np.exp((1 - self.alpha) * self.F * eta / (self.R * self.T))
        term_reverse = np.exp(-self.alpha * self.F * eta / (self.R * self.T))
        return self.j0 * (term_forward - term_reverse)

    def compute_butler_volmer_curve(self, 
                                    eta_min=Config.BV_ETA_MIN,
                                    eta_max=Config.BV_ETA_MAX,
                                    num_points=Config.BV_NUM_PTS
                                    ):
        eta_range = np.linspace(eta_min, eta_max, num_points)
        j_values = self.butler_volmer(eta_range)
        return eta_range, j_values

    ###########################################################################
    # 5. Marcus Theory (All in eV)
    ###########################################################################
    def marcus_rate(self, deltaG_eV: float) -> float:
        denom = 4.0 * max(self.lambda_marcus, Config.EPSILON_CONST) * Config.BOLTZMANN_EVK_CONST * self.T + Config.EPSILON_CONST
        numerator = (self.lambda_marcus + deltaG_eV)**2
        exponent = - (numerator / denom)
        return self.A_marcus * np.exp(exponent)

    def compute_marcus_curve(self, 
                            deltaG0_eV: float = None,
                            window_eV=Config.MARCUS_WIN_EV,
                            num_points=Config.MARCUS_NUM_PTS
                            ):
        if deltaG0_eV is None:
            deltaG0_eV = self.deltaG0
        deltaG_min = deltaG0_eV - window_eV
        deltaG_max = deltaG0_eV + window_eV
        deltaG_range = np.linspace(deltaG_min, deltaG_max, num_points)
        k_et = np.array([self.marcus_rate(dG) for dG in deltaG_range])
        return deltaG_range, k_et

    ###########################################################################
    # 6. Integrated Model
    ###########################################################################
    def compute_integrated_current_curve(self, 
                                        method='butler-volmer', 
                                        t_min=Config.BV_T_MIN,
                                        t_max=Config.BV_T_MAX,
                                        num_points=Config.BV_NUM_PTS
                                        ):
        Kd = self.k_off / (self.k_on + Config.EPSILON_CONST)
        T_conc_range = np.logspace(np.log10(t_min), np.log10(t_max), num_points)
        fraction_bound = self.hill_langmuir(T_conc_range, Kd, self.n)
        integrated_current = np.zeros_like(T_conc_range)
        
        if method == 'butler-volmer':
            # reference
            j_max = self.butler_volmer(self.eta0)
            for i, theta in enumerate(fraction_bound):
                eta_eff = self.eta0 * (1 - theta)
                j = self.butler_volmer(eta_eff)
                j_norm = j / (j_max + Config.EPSILON_CONST)
                integrated_current[i] = self.I_min + (self.I_max - self.I_min) * j_norm

        elif method == 'marcus':
            # reference
            k_max = self.marcus_rate(self.deltaG0)
            for i, theta in enumerate(fraction_bound):
                effective_deltaG = self.deltaG0 - self.deltaG_shift * theta
                k_et = self.marcus_rate(effective_deltaG)
                k_norm = k_et / (k_max + Config.EPSILON_CONST)
                integrated_current[i] = self.I_min + (self.I_max - self.I_min) * k_norm

        else:
            # fallback: linear HL response
            integrated_current = self.I_min + (self.I_max - self.I_min) * fraction_bound
        
        return T_conc_range, fraction_bound, integrated_current
