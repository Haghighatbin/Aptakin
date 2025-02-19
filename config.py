class Config:
    EPSILON_CONST       = 1e-12             # To avoid division by zero error
    FARADAY_CONST       = 96485.3329        # Faraday's constant (C/mol)
    UNI_GAS_CONST       = 8.314             # Universal gas constant (J/(mol·K))
    BOLTZMANN_JK_CONST  = 1.380649e-23      # Boltzmann constant (J/K)
    BOLTZMANN_EVK_CONST = 8.617333262145e-5 # Boltzmann constant (eV/K)

    # -------------------------------------------------
    # 1. Binding Kinetics Parameters
    # -------------------------------------------------
    BIND_KIN_APT_CONC_TOTAL = 1.0           # Total aptamer concentration (arbitrary units)
    BIND_KIN_TRGT_CONC      = 10.0          # Target concentration
    BIND_KIN_K_ON           = 1.0           # Binding rate constant (1/s)
    BIND_KIN_K_OFF          = 0.5           # Unbinding rate constant (1/s)
    BIND_KIN_I_MIN          = 0.1           # Minimum current
    BIND_KIN_I_MAX          = 1.0           # Maximum current
    BIND_KIN_APT_INIT_CONC  = 0.0           # Initial bound aptamer concentration
    BIND_KIN_HILL_COEFF     = 1.0           # Hill coefficient
    BIND_KIN_T_END          = 4
    BIND_KIN_NUM_PTS        = 200

    # -------------------------------------------------
    # 2. SWV Simulation Parameters
    # -------------------------------------------------
    SWV_SIM_K_ET_UNBOUND    = 10
    SWV_SIM_K_ET_BOUND      = 50
    SWV_SIM_FREQ_DECAY      = 100
    SWV_SIM_I0              = 1.0
    SWV_PARAM_NUM_PTS       = 200

    # -------------------------------------------------
    # 3. Butler–Volmer Parameters
    # -------------------------------------------------
    BV_CURR_DENS_J0         = 1e-3          # Exchange current density (A/cm²)
    BV_ALPHA                = 0.5           # Charge transfer coefficient
    BV_TEMP                 = 298.15        # Temperature (K)
    BV_ETA_0                = 0.5
    BV_ETA_MIN              = -0.5
    BV_ETA_MAX              = 0.5
    BV_NUM_PTS              = 200
    BV_T_MIN                = 0.01
    BV_T_MAX                = 10
    BV_NUM_PTS              = 1000  

    # -------------------------------------------------
    # 4. Marcus Theory Parameters (energy in eV)
    # -------------------------------------------------
    MARCUS_PREEXP_A         = 1e12          # Pre-exponential factor (s⁻¹)
    MARCUS_REORG_E_LAMBDA   = 0.2           # Reorganisation energy (eV)
    MARCUS_DELTAG           = -0.1          # Standard free energy change (eV)
    MARCUS_DELTAG_SHIFT     = 0.1
    MARCUS_WIN_EV           = 0.5
    MARCUS_NUM_PTS          = 200

    # -------------------------------------------------
    # 5.Experimental Simulation Paramters
    # -------------------------------------------------
    EXPRM_DATA_FREQ_LOW          = 16      # Frequency below the critical frequency (e.g. 16 Hz)
    EXPRM_DATA_FREQ_HIGH         = 240     # Frequency above the critical frequency (e.g. 240 Hz)
    EXPRM_DATA_CHRG_LOW_UNBOUND  = 5.0
    EXPRM_DATA_CHRG_LOW_BOUND    = 5.0
    EXPRM_DATA_CHRG_HIGH_UNBOUND = 2.0
    EXPRM_DATA_CHRG_HIGH_BOUND   = 2.0
    EXPRM_DATA_ELEC_AREA         = 15.0    # Electrode geometry (nominal area in mm²)
    EXPRM_DATA_MB_N_ET           = 2       # Number of electrons transferred during MB redox (default 2)

    # -------------------------------------------------
    # 6.Hill-Langmuir Isotherm Paramters
    # -------------------------------------------------
    HL_ISO_T_MIN                 = 0.01
    HL_ISO_T_MAX                 = 10
    HL_ISO_NUM_PTS               = 1000 

    # -------------------------------------------------
    # Layout and Callback Settings
    # -------------------------------------------------
    PLOT_TITLE_FONTSIZE = 20
    PLOT_XAXIS_FONTSIZE = 20
    PLOT_YAXIS_FONTSIZE = 20 
    PLOT_TICK_FONTSIZE  = 18