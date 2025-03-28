class Config:
    # 1. BOUND-STATE MODEL PARAMETERS
    BOUND_K_ET = 20.0                                           # Electron-transfer rate proxy in exponential model.
    BOUND_F_DECAY1 = 20.0                                       # First characteristic frequency decay.
    BOUND_F_DECAY2 = 50.0                                       # Second characteristic frequency decay.
    BOUND_A1 = 0.7                                              # Weight factor for double-exponential (if parallel sum is used).

    # 2. UNBOUND-STATE MODEL PARAMETERS
    UNBOUND_K_ET = 20                                           # Electron-transfer rate proxy in exponential model.
    UNBOUND_F_DECAY1 = 20.0                                     # First characteristic frequency decay.
    UNBOUND_F_DECAY2 = 50.0                                     # Second characteristic frequency decay.
    UNBOUND_A1 = 0.7                                            # Weight factor for double-exponential (if parallel sum is used).

    # 3. CURVE-FITTING SETTINGS FOR THE BOUND STATE
    CF_BOUND_MAXFEV = 10000                                     # Maximum function evaluations for curve_fit.
    CF_BOUND_I0_MUL = 10                                        # Multiplier for initial guess of I0 (bound).
    CF_BOUND_L_BNDS = [0,    0,    0,    0]                     # Lower bounds for [I0, k_et, f_decay1, f_decay2].
    CF_BOUND_H_BNDS = [1e-4, 1e3, 1e3, 1e3]                     # Upper bounds for [I0, k_et, f_decay1, f_decay2].
    CF_BOUND_PCT_SIGMA = 0.1                                    # Sigma settings for weighting data points in curve_fit: 0.1 means 10% relative error.
    CF_BOUND_RAW_VAL_SIGMA = 5e-8                               # Absolute error fallback if needed.
    CF_BOUND_METHOD = 'trf'                                     # Method for curve_fit; "trf" often used for bounded problems.

    # 4. CURVE-FITTING SETTINGS FOR THE UNBOUND STATE
    CF_UNBOUND_MAXFEV = 10000                                   # Maximum function evaluations for curve_fit.
    CF_UNBOUND_I0_MUL = 1                                       # Multiplier for initial guess of I0 (unbound).
    CF_UNBOUND_L_BNDS = [0,    0,    0,    0]                   # Lower bounds for [I0, k_et, f_decay1, f_decay2].
    CF_UNBOUND_H_BNDS = [1e-3, 1e4, 1e4, 1e4]                   # Upper bounds for [I0, k_et, f_decay1, f_decay2].
    CF_UNBOUND_PCT_SIGMA = 0.1                                  # Sigma settings for weighting data points in curve_fit: 0.1 means 10% relative error.
    CF_UNBOUND_RAW_VAL_SIGMA = 5e-8                             # Absolute error fallback if needed.
    CF_UNBOUND_METHOD = 'trf'                                   # Method for curve_fit; "trf" often used for bounded problems.

    # 5. PLOTTING PARAMETERS (EXPERIMENTAL DATA & EXTRAPOLATION)
    PLT_LNSPC_NUM = 200                                         # Number of points in linspace for plotting fits.
    PLT_NODEVZERO_VAL = 1e-12                                   # Small constant to avoid division by zero.
    PLT_EXTP_MAXFREQ = 600                                      # Upper frequency bound for extrapolated plotting.
    PLT_EXTP_MINFREQ = 2                                        # Lower frequency bound for extrapolated plotting.

    # 6. PHYSICAL CONSTANTS
    FARADAY_CONST       = 96485.3329                            # Faraday's constant (C/mol).
    UNI_GAS_CONST       = 8.314                                 # Universal gas constant (J/(mol·K)).
    BOLTZMANN_JK_CONST  = 1.380649e-23                          # Boltzmann constant (J/K).
    BOLTZMANN_EVK_CONST = 8.617333262145e-5                     # Boltzmann constant (eV/K).
    STRD_TEMP           = 298.15                                # Standard temperature in Kelvin (25 °C).

    # 7. BUTLER–VOLMER PARAMETERS
    BV_CURR_DENS_J0     = 1e-3                                  # Exchange current density j0 (A/cm²).
    BV_ALPHA            = 0.5                                   # Charge transfer coefficient (typical range 0.2–0.6).
    BV_ETA_0            = 0.5                                   # Nominal overpotential in Volts (if used as constant).
    BV_ETA_MIN          = -0.5                                  # Lower bound for scanning overpotential (if needed).
    BV_ETA_MAX          = 0.5                                   # Upper bound for scanning overpotential.
    BV_NUM_PTS          = 1000                                  # Number of points for plotting BV curves.
    BV_T_MIN            = 0.01                                  # 
    BV_T_MAX            = 10                                    # 

    # Frequency-dependent "effective" overpotential factor:
    BV_UNBOUND_K_EFF    = 5                                     # Rate at which overpotential transitions for unbound.
    BV_BOUND_K_EFF      = 3                                     # Rate at which overpotential transitions for bound.

    BV_UNBOUND_F_DECAY1 = 20                                    # Decay param #1 for unbound (double-exponential).
    BV_UNBOUND_F_DECAY2 = 50                                    # Decay param #2 for unbound.
    BV_BOUND_F_DECAY1   = 20                                    # Decay param #1 for bound.
    BV_BOUND_F_DECAY2   = 50                                    # Decay param #2 for bound.

    # 8. MARCUS THEORY PARAMETERS
    MARCUS_PREEXP_A         = 1e12                              # Pre-exponential factor (s⁻¹) in Marcus ET rate.
    MARCUS_REORG_E_LAMBDA   = 0.2                               # Reorganisation energy λ (eV).
    MARCUS_DELTAG           = -0.1                              # Standard free energy change ΔG⁰ (eV).
    MARCUS_DELTAG_SHIFT     = 0.1                               # Additional shift if needed.
    MARCUS_WIN_EV           = 0.5                               # Range (±) around ΔG⁰ for plotting Marcus curves.
    MARCUS_NUM_PTS          = 200                               # Points for Marcus curve plots.

    # Frequency-dependent factor for Marcus-based SWV:
    MARCUS_UNBOUND_K_FREQ   = 5                                 # Controls how quickly unbound ET rate saturates with freq.
    MARCUS_BOUND_K_FREQ     = 5                                 # Same but for bound state.
    
    MARCUS_UNBOUND_F_DECAY1 = 20                                # Decay param #1 (unbound).
    MARCUS_UNBOUND_F_DECAY2 = 50                                # Decay param #2 (unbound).
    MARCUS_BOUND_F_DECAY1   = 20                                # Decay param #1 (bound).
    MARCUS_BOUND_F_DECAY2   = 50                                # Decay param #2 (bound).
