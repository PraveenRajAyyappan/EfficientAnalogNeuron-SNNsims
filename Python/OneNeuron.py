# Written by Jaron Rosenberg
# This script simulates the response of a singular FH type neuron and its synapse
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from NeuronModules import *

# Define Simulation Params and Run
dt = 0.00001  # In s
tspan = [0, 0.2]  # In s
Params['tvec'] = np.arange(tspan[0], tspan[1] + dt, dt)

# Run a quick transient sim with zero input for the true initial conditions
y0 = np.array(Params['DC_OFF'])
Params['Input'] = np.zeros(len(Params['tvec']))  # In pA
y = odeint(NeuronTest, y0, Params['tvec'], args=(Params,))

# Run the actual sim
Params['Input'] = 35 * ((Params['tvec'] >= 0.050) & (Params['tvec'] <= 0.15))  # In pA
y = odeint(NeuronTest, y[-1, :], Params['tvec'], args=(Params,))

# Plot Na Clamp Data
# Plot
fig, ax = plt.subplots(5, 1, figsize=(6, 10), sharex=True)
ax[0].plot(Params['tvec'] * 1000, y[:, 0], 'k')
ax[0].set_ylabel(r'$V_{mem}$ $(U_T)$')
ax[1].plot(Params['tvec'] * 1000, y[:, 1], 'k')
ax[1].set_ylabel(r'$V_{K}$ $(U_T)$')
ax[2].plot(Params['tvec'] * 1000, y[:, 2], 'k')
ax[2].set_ylabel(r'$V_{g}$ $(U_T)$')
ax[3].plot(Params['tvec'] * 1000, y[:, 3], 'k')
ax[3].set_ylabel(r'$V_{Na}$ $(U_T)$')
ax[4].plot(Params['tvec'] * 1000, Params['Input'], 'k')
ax[4].set_ylabel(r'$I_{in}$ $(pA)$')
ax[4].set_xlabel('Time (ms)')

plt.tight_layout()
plt.xlim([tspan[0] * 1000, tspan[1] * 1000])
plt.show()
