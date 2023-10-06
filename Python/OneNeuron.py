# Written by Jaron Rosenberg
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint

# Neuron Computation
def Neuron(y, t, Params):
    # FH Type HH Neuron ODE
    Kp = 0.75
    I_Vec = [
        np.interp(t, Params['tvec'], Params['Input']),
        1 - Params['EXPMEL'] * np.exp(y[0]),
        np.exp(-Kp * y[3]) * (1 - np.exp(y[0] - Params['ENA'])),
        np.exp(-Kp * y[2]),
        1 - Params['EXPVSAT'] * np.exp(-y[3]),
        np.exp(y[0] - Kp * y[1]) - np.exp(Params['EK'] - Kp * y[1]),
        1 - Params['EXPMVGK'] * np.exp(y[1]),
        np.exp(y[2]) - np.exp(y[3])
    ]
    
    dydt = np.dot(Params['Neuron_Mat'], I_Vec)
    return dydt


# Initialize
Params = {}

# Define Neuron (Physical) Params
# Neuron
Params['CMEM'] = 800e-15
Params['CZ'] = 120e-15
Params['CNA'] = 960e-15
Params['CK'] = 400e-15
Params['CW'] = 8e-15
Params['CL'] = 8e-15
Params['M_Na'] = 1 / 0.6  # Microns
Params['M_K'] = 1 / 0.6  # Microns
Params['M_Tn'] = 1 / 0.6  # Microns
Params['M_Tm'] = 12 / 0.6  # Microns
Params['M_Amp'] = 12 / 1  # Microns
Params['M_Th'] = 1 / 0.6  # Microns
Params['M_Tl'] = 1 / 0.6  # Microns

# Define Voltage Params
# Misc
UT = 25.0e-3  # Thermal voltage
Params['DC_OFF'] = [46, 14, 26.55, 26.55]  # Estimated initial DC Offsets (in UT). Does not need to be perfect.

# Neuron
Params['ENA'] = 1.30 / UT
Params['VAMP'] = 1.440 / UT
Params['VTH'] = -0.27 / UT
Params['VTM'] = 0.7125 / UT
Params['VSAT'] = 0.385 / UT
Params['VTN'] = -0.68 / UT
Params['VGK'] = 0.350 / UT
Params['EK'] = 1.15 / UT
Params['VTL'] = 0.45 / UT
Params['EL'] = 1.15 / UT

# Transfer Matrices
# Transistor Parameters
VDD = 100
Kp = 0.75
Ithp = 125e-9
VT0p = 28
Kn = 0.67
Ithn = 300e-9
VT0n = 20

# Define capacitance elements
C_A_Square = (Params['CNA'] + Params['CZ'] + Params['CW']) * (Params['CZ'] + Params['CL']) - Params['CZ']**2
C_B_Cube = (Params['CMEM'] + Params['CNA']) * C_A_Square - (Params['CNA']**2) * (Params['CZ'] + Params['CL'])
C_11 = C_A_Square / C_B_Cube
C_14 = Params['CNA'] * Params['CZ'] / C_B_Cube
C_18 = Params['CNA'] * Params['CL'] / C_B_Cube
C_27 = C_A_Square / C_B_Cube + 1 / Params['CK']
C_31 = Params['CNA'] * (Params['CZ'] + Params['CL']) / C_B_Cube
C_34 = Params['CZ'] / C_A_Square + (Params['CNA']**2) * Params['CZ'] * (Params['CZ'] + Params['CL']) / (C_A_Square * C_B_Cube)
C_38 = Params['CL'] / C_A_Square + (Params['CNA']**2) * Params['CL'] * (Params['CZ'] + Params['CL']) / (C_A_Square * C_B_Cube)
C_44 = (Params['CNA']**2 * Params['CZ']**2) / (C_A_Square * C_B_Cube) + (Params['CZ']**2) / (C_A_Square * (Params['CZ'] + Params['CL'])) + 1 / (Params['CZ'] + Params['CL'])
C_48 = (Params['CL'] * Params['CZ'] - C_A_Square) / (C_A_Square * (Params['CZ'] + Params['CL'])) + (Params['CNA']**2 * Params['CL'] * Params['CZ']) / (C_A_Square * C_B_Cube)

# Compute Neuron Transfer Matrix
Neuron_Mat = (1 / UT) * np.array([[C_11, C_11, C_11, C_14, -C_14, -C_11, -C_11, -C_18],
                                  [C_11, C_11, C_11, C_14, -C_14, -C_11, -C_27, -C_18],
                                  [C_31, C_31, C_31, C_34, -C_34, -C_31, -C_31, C_38],
                                  [C_14, C_14, C_14, C_44, -C_44, -C_14, -C_14, -C_48]])

repmat_vector = np.array([1e-12, Params['M_Tl'] * Ithp * np.exp(Kp * (VDD - Params['VTL'] - VT0p) - VDD + Params['EL']),
                           Params['M_Na'] * Ithp * np.exp(Kp * (VDD - VT0p) - VDD + Params['ENA']),
                           Params['M_Amp'] * Ithp * np.exp(Kp * (VDD - VT0p) - VDD + Params['VAMP']),
                           Params['M_Tm'] * Ithn * np.exp(Kn * (Params['VTM'] - VT0n) - Params['VSAT']),
                           Params['M_K'] * Ithp * np.exp(Kp * (VDD - VT0p) - VDD),
                           -Params['M_Tn'] * Ithp * np.exp(Kp * (VDD - VT0p - Params['VTN']) - VDD + Params['VGK']),
                           Params['M_Th'] * Ithp * np.exp(Kp * (VDD - Params['VTH'] - VT0p) - VDD)])

Params['Neuron_Mat'] = Neuron_Mat * np.tile(repmat_vector, (4, 1))
Params['EXPMEL'] = np.exp(-Params['EL'])
Params['EXPMVGK'] = np.exp(-Params['VGK'])
Params['EXPVSAT'] = np.exp(Params['VSAT'])

# Define Simulation Params and Run
dt = 0.00001  # In s
tspan = [0, 0.2]  # In s
Params['tvec'] = np.arange(tspan[0], tspan[1] + dt, dt)

# Run a quick transient sim with zero input for the true initial conditions
y0 = np.array(Params['DC_OFF'])
Params['Input'] = np.zeros(len(Params['tvec']))  # In pA
y = odeint(Neuron, y0, Params['tvec'], args=(Params,))

# Run the actual sim
Params['Input'] = 35 * ((Params['tvec'] >= 0.050) & (Params['tvec'] <= 0.15))  # In pA
y = odeint(Neuron, y[-1, :], Params['tvec'], args=(Params,))

# Plot Na Clamp Data
# Plot
fig, ax = plt.subplots(5, 1, figsize=(6, 10), sharex=True)
ax[0].plot(Params['tvec'] * 1000, y[:, 0], 'k')
ax[0].set_ylabel('V_{mem} (U_T)')
ax[1].plot(Params['tvec'] * 1000, y[:, 1], 'k')
ax[1].set_ylabel('V_{K} (U_T)')
ax[2].plot(Params['tvec'] * 1000, y[:, 2], 'k')
ax[2].set_ylabel('V_{g} (U_T)')
ax[3].plot(Params['tvec'] * 1000, y[:, 3], 'k')
ax[3].set_ylabel('V_{Na} (U_T)')
ax[4].plot(Params['tvec'] * 1000, Params['Input'], 'k')
ax[4].set_ylabel('I_{in} (pA)')
ax[4].set_xlabel('Time (ms)')

plt.tight_layout()
plt.xlim([tspan[0] * 1000, tspan[1] * 1000])
plt.show()
