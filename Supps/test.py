import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint

# Define the model parameters
# A_total = 1.0      # Total aptamer concentration (normalized)
# T_conc = 10.0      # Target concentration (assumed constant and in excess)
# k_on = 1.0         # Binding rate constant
# k_off = 0.5        # Unbinding rate constant

# # Sensor current parameters (arbitrary units)
# I_min = 0.1        # Current when no aptamers are bound
# I_max = 1.0        # Current when all aptamers are bound

# # Initial condition: at time t=0, no aptamer is bound
# AT0 = 0.0

# # Define the differential equation for [AT]
# def dAT_dt(AT, t, k_on, k_off, A_total, T_conc):
#     """
#     Compute the time derivative of the bound aptamer concentration.
    
#     Parameters:
#     - AT: current bound aptamer concentration
#     - t: time (unused explicitly as the equation is time-invariant)
#     - k_on, k_off: kinetic rate constants
#     - A_total: total aptamer concentration
#     - T_conc: target concentration (assumed constant)
    
#     Returns:
#     - dAT/dt: time derivative of bound aptamer concentration
#     """
#     A_free = A_total - AT  # Aptamer that is not bound
#     dAT = k_on * A_free * T_conc - k_off * AT
#     return dAT

# # Set up the time vector for the simulation
# t = np.linspace(0, 10, 200)  # Simulate from 0 to 10 seconds

# # Solve the differential equation using odeint
# AT = odeint(dAT_dt, AT0, t, args=(k_on, k_off, A_total, T_conc))
# AT = AT.flatten()  # Flatten the result for easier processing

# # Calculate the fraction of bound aptamers
# fraction_bound = AT / A_total

# # Calculate the sensor current based on the fraction bound
# I = I_min + (I_max - I_min) * fraction_bound

# Plot the fraction bound over time
# plt.figure(figsize=(8, 4))
# plt.plot(t, fraction_bound, label='Fraction Bound', color='blue')
# plt.xlabel('Time (s)')
# plt.ylabel('Fraction Bound')
# plt.title('Aptamer Binding Kinetics')
# plt.legend()
# plt.grid(True)
# plt.show()

# # Plot the sensor current response over time
# plt.figure(figsize=(8, 4))
# plt.plot(t, I, label='Sensor Current', color='red')
# plt.xlabel('Time (s)')
# plt.ylabel('Current (a.u.)')
# plt.title('Electrochemical Sensor Response')
# plt.legend()
# plt.grid(True)
# plt.show()

def hill_langmuir(T, Kd, n):
    """
    Compute the fraction of bound aptamer using the Hill-Langmuir isotherm.
    
    Parameters:
    - T: target concentration (can be a scalar or numpy array)
    - Kd: dissociation constant
    - n: Hill coefficient
    
    Returns:
    - fraction_bound: Fraction of aptamers bound
    """
    return T**n / (Kd**n + T**n)

# Example parameters
# Kd = k_off / k_on  # Dissociation constant from our kinetic model
# n = 1.0           # Hill coefficient (set to 1 for non-cooperative binding)
# T_conc_range = np.linspace(0, 20, 200)  # Range of target concentrations

# # Compute fraction bound using Hill-Langmuir isotherm
# fraction_bound_hl = hill_langmuir(T_conc_range, Kd, n)

# # Compute sensor current based on fraction bound
# I_hl = I_min + (I_max - I_min) * fraction_bound_hl

# # Plot the Hill-Langmuir binding curve
# plt.figure(figsize=(8, 4))
# plt.plot(T_conc_range, fraction_bound_hl, label='Hill-Langmuir Fit (n=1)', color='green')
# plt.xlabel('Target Concentration (T)')
# plt.ylabel('Fraction Bound')
# plt.title('Equilibrium Binding Curve')
# plt.legend()
# plt.grid(True)
# plt.show()

# # Plot the sensor current response based on the Hill-Langmuir fit
# plt.figure(figsize=(8, 4))
# plt.plot(T_conc_range, I_hl, label='Sensor Current (Hill-Langmuir)', color='magenta')
# plt.xlabel('Target Concentration (T)')
# plt.ylabel('Current (a.u.)')
# plt.title('Sensor Response Curve')
# plt.legend()
# plt.grid(True)
# plt.show()


def hill_langmuir(T, Kd, n):
    """
    Compute the fraction of bound aptamer using the Hill-Langmuir isotherm.
    
    Parameters:
    - T: target concentration (scalar or numpy array)
    - Kd: dissociation constant
    - n: Hill coefficient
    
    Returns:
    - Fraction bound: [T]^n / (Kd^n + [T]^n)
    """
    return T**n / (Kd**n + T**n)

# Parameters
# Kd = 5.0  # Example dissociation constant
# T = np.linspace(0.1, 50, 300)  # Target concentration range

# # Case 1: n = 1 (non-cooperative binding)
# fraction_bound_n1 = hill_langmuir(T, Kd, 1)

# # Case 2: n = 2 (cooperative binding, resulting in a more sigmoidal shape)
# fraction_bound_n2 = hill_langmuir(T, Kd, 2)

# # Plot on a linear concentration axis
# plt.figure(figsize=(8, 4))
# plt.plot(T, fraction_bound_n1, label='n = 1 (Hyperbolic)', color='blue')
# plt.plot(T, fraction_bound_n2, label='n = 2 (Sigmoidal)', color='red')
# plt.xlabel('Target Concentration [T]')
# plt.ylabel('Fraction Bound')
# plt.title('Hill-Langmuir Isotherm (Linear Scale)')
# plt.legend()
# plt.grid(True)
# plt.show()

# Plot on a logarithmic concentration axis
# plt.figure(figsize=(8, 4))
# plt.semilogx(T, fraction_bound_n1, label='n = 1 (Appears Sigmoidal)', color='blue')
# plt.semilogx(T, fraction_bound_n2, label='n = 2 (Sigmoidal)', color='red')
# plt.xlabel('Target Concentration [T] (log scale)')
# plt.ylabel('Fraction Bound')
# plt.title('Hill-Langmuir Isotherm (Logarithmic Scale)')
# plt.legend()
# plt.grid(True, which="both", ls="--")
# plt.show()


def swv_pulse_current(I0, k_et, pulse_duration):
    """
    Simulate the current at the end of a pulse,
    assuming an exponential decay in response to a potential step.
    
    Parameters:
      I0            : Initial current amplitude at the start of the pulse.
      k_et          : Electron-transfer rate constant (s^-1).
      pulse_duration: Duration of the pulse (s).
      
    Returns:
      I_pulse: Current at the end of the pulse.
    """
    return I0 * np.exp(-k_et * pulse_duration)

def simulate_swv(frequency, k_et, I0=1.0):
    """
    Simulate a simplified net SWV current for one frequency.
    
    Parameters:
      frequency: Frequency of the SWV (Hz).
      k_et     : Electron-transfer rate constant (s^-1) for the state (bound or unbound).
      I0       : Initial current (arbitrary units).
    
    Returns:
      I_net: Simulated net SWV current for that frequency.
    """
    # Calculate the period T (in seconds) from the frequency
    T = 1 / frequency
    
    # Each square-wave pulse lasts half of the period:
    pulse_duration = T / 2.0  # Duration of each forward or reverse pulse
    
    # Simulate the current at the end of the forward pulse
    I_forward = swv_pulse_current(I0, k_et, pulse_duration)
    
    # Simulate the current at the end of the reverse pulse.
    # In a real SWV, these pulses are at different potentials, so the currents differ.
    # Here, we introduce a small asymmetry (multiply by 0.95) for illustration.
    I_reverse = swv_pulse_current(I0, k_et, pulse_duration)
    
    # The net SWV current is the difference between the forward and reverse responses.
    I_net = I_forward - I_reverse * 0.95
    
    return I_net

# # Define a range of frequencies to simulate (in Hz)
# frequencies = np.linspace(1, 1000, 200)

# # Define electron-transfer rate constants for the two states:
# # Unbound state: fast electron transfer (higher k_et)
# k_et_unbound = 50   # s^-1
# # Bound state: slower electron transfer due to conformational change (lower k_et)
# k_et_bound = 10     # s^-1

# # Simulate the net SWV current over the frequency range for both states.
# I_net_unbound = [simulate_swv(f, k_et_unbound) for f in frequencies]
# I_net_bound   = [simulate_swv(f, k_et_bound)   for f in frequencies]

# # Plot the simulated SWV responses versus frequency.
# plt.figure(figsize=(8, 5))
# plt.plot(frequencies, I_net_unbound, label='Unbound (fast k_et)', color='blue')
# plt.plot(frequencies, I_net_bound, label='Bound (slow k_et)', color='red')
# plt.xlabel('Frequency (Hz)')
# plt.ylabel('Simulated Net SWV Current (a.u.)')
# plt.title('Simulated SWV Response vs Frequency')
# plt.legend()
# plt.grid(True)
# plt.show()


# def swv_pulse_current(I0, k_et, pulse_duration):
#     """
#     Compute the current at the end of a square-wave pulse,
#     assuming an exponential decay due to electron-transfer kinetics.
    
#     Parameters:
#       I0            : Initial current at the beginning of the pulse.
#       k_et          : Electron-transfer rate constant (s^-1).
#       pulse_duration: Duration of the pulse (s).
      
#     Returns:
#       I_pulse: Current at the end of the pulse.
#     """
#     return I0 * np.exp(-k_et * pulse_duration)

# def simulate_swv(frequency, k_et, I0=1.0, f_decay=100):
#     """
#     Simulate the net SWV current at a given frequency using a model that 
#     includes:
#       - An exponential term for the decay due to electron-transfer kinetics.
#       - An additional exponential term that causes a drop at high frequency.
    
#     Parameters:
#       frequency: SWV frequency in Hz.
#       k_et: Electron-transfer rate constant (s^-1).
#       I0: Initial current amplitude (arbitrary units).
#       f_decay: Frequency decay constant (Hz), which controls the high-frequency drop.
    
#     Returns:
#       I_net: Net SWV current (arbitrary units).
#     """
#     return I0 * np.exp(-k_et/(2*frequency)) * np.exp(-frequency/f_decay)

# # Define frequency range for the simulation (Hz)
# frequencies = np.linspace(1, 1000, 200)

# # Define parameters for the two states:
# # Unbound state: slower kinetics (MB is farther from the electrode)
# k_et_unbound = 8   # s^-1
# # Bound state: faster kinetics (MB is closer to the electrode)
# k_et_bound   = 80   # s^-1
# I0 = 1.0
# f_decay = 100  # Hz

# # Compute the net SWV current over the frequency range for each state.
# I_net_unbound = np.array([simulate_swv(f, k_et_unbound, I0, f_decay) for f in frequencies])
# I_net_bound   = np.array([simulate_swv(f, k_et_bound, I0, f_decay) for f in frequencies])

# # Determine the critical frequencies (frequency at maximum net current)
# critical_freq_unbound = frequencies[np.argmax(I_net_unbound)]
# critical_freq_bound   = frequencies[np.argmax(I_net_bound)]

# print(f"Critical frequency (unbound): {critical_freq_unbound:.1f} Hz")
# print(f"Critical frequency (bound):   {critical_freq_bound:.1f} Hz")

# # Plot the simulated SWV responses versus frequency.
# plt.figure(figsize=(10, 6))
# plt.plot(frequencies, I_net_unbound, label=f'Unbound (k_et={k_et_unbound} s⁻¹)', color='blue')
# plt.plot(frequencies, I_net_bound, label=f'Bound (k_et={k_et_bound} s⁻¹)', color='red')
# plt.axvline(critical_freq_unbound, color='blue', linestyle='--',
#             label=f'Critical f (unbound) = {critical_freq_unbound:.1f} Hz')
# plt.axvline(critical_freq_bound, color='red', linestyle='--',
#             label=f'Critical f (bound) = {critical_freq_bound:.1f} Hz')
# plt.xlabel('Frequency (Hz)')
# plt.ylabel('Simulated Net SWV Current (a.u.)')
# plt.title('Simulated SWV Response vs Frequency with Critical Frequencies')
# plt.legend()
# plt.grid(True)
# plt.show()

import numpy as np
import matplotlib.pyplot as plt

def simulate_swv(frequency, k_et, I0=1.0, f_decay=100):
    """
    Simulate the net SWV current at a given frequency using a model that 
    includes:
      - An exponential term for the decay due to electron-transfer kinetics.
      - An additional exponential term that causes a drop at high frequency.
    
    Parameters:
      frequency: SWV frequency in Hz.
      k_et: Electron-transfer rate constant (s^-1).
      I0: Initial current amplitude (arbitrary units).
      f_decay: Frequency decay constant (Hz), which controls the high-frequency drop.
    
    Returns:
      I_net: Net SWV current (arbitrary units).
    """
    return I0 * np.exp(-k_et/(2*frequency)) * np.exp(-frequency/f_decay)

# Define frequency range for the simulation (Hz)
frequencies = np.linspace(1, 1000, 200)

# Define parameters for the two states:
# Unbound state: slower kinetics (MB is farther from the electrode)
k_et_unbound = 10   # s^-1
# Bound state: faster kinetics (MB is closer to the electrode)
k_et_bound   = 50   # s^-1
I0 = 1.0
f_decay = 100  # Hz

# Compute the net SWV current over the frequency range for each state.
I_net_unbound = np.array([simulate_swv(f, k_et_unbound, I0, f_decay) for f in frequencies])
I_net_bound   = np.array([simulate_swv(f, k_et_bound, I0, f_decay) for f in frequencies])

# Normalize by their maximum values to compare their relative shapes
I_net_unbound /= np.max(I_net_unbound)
I_net_bound   /= np.max(I_net_bound)

# Find the crossover frequency (where the two normalized curves intersect)
crossover_index = np.argmin(np.abs(I_net_unbound - I_net_bound))
crossover_freq = frequencies[crossover_index]

print(f"Crossover Frequency: {crossover_freq:.1f} Hz")

# Plot the simulated SWV responses versus frequency.
plt.figure(figsize=(10, 6))
plt.plot(frequencies, I_net_unbound, label=f'Unbound (k_et={k_et_unbound} s⁻¹) (Normalized)', color='blue')
plt.plot(frequencies, I_net_bound, label=f'Bound (k_et={k_et_bound} s⁻¹) (Normalized)', color='red')
plt.axvline(crossover_freq, color='black', linestyle='--',
            label=f'Crossover f = {crossover_freq:.1f} Hz')
plt.xlabel('Frequency (Hz)')
plt.ylabel('Normalized Net SWV Current')
plt.title('Normalized SWV Response vs Frequency with Crossover Point')
plt.legend()
plt.grid(True)
plt.show()
