# This Script Simulates the response of a Farquhar-Hasler Neuron in a WTA
import numpy as np
import matplotlib.pyplot as plt
from NeuronModules import *


# Define Network Connectivity Via the Synapse Weight Matrix
# Note: Odd number columns are excitatory, even number columns are inhibitory
Params['NeuronPopulation'] = 4  # Number of neurons in the network
Ex = 0.26  # Excitatory synapse "strength" (Lower value means greater weightage)
Ih = 0.20  # Inhibitory synapse "strength" (Lower value means greater weightage)
Weights = np.array([[5, 5, 5, 5, 5, 5, 5, Ih],
                    [5, 5, 5, 5, 5, 5, 5, Ih],
                    [5, 5, 5, 5, 5, 5, 5, Ih],
                    [Ex, 5, Ex, 5, Ex, 5, 5, 5]])  # Controls neuron connectivity
Params['ExternalCurrentIdx'] = [1, 2, 3]  # External current will be channeled into neurons with these indices
Params['Synapse_Mat'] = Params['M_Syn'] * Ithp * np.exp(Kp * (VDD - VT0p) - VDD) * \
    np.exp(-Kp * Weights / UT) * ((Weights < 2.5))  # Computes the final transfer matrix

# Estimate y0 by producing an initial guess
y0 = np.array([])
for i in range(Params['NeuronPopulation']):
    y0 = np.append(y0, np.array([100, 100]))
for i in range(Params['NeuronPopulation']):
    y0 = np.append(y0, np.array(Params['DC_OFF']))

# Define Simulation Params and Run
dt =  0.00001  # In s
tspan = [0, 0.8]  # In s
Params['tvec'] = np.arange(tspan[0], tspan[1] + dt, dt)

# Run a quick transient sim with zero input for the true initial conditions
Params['Input'] = np.zeros((len(Params['tvec']), 3))  # In pA
y = euler_solver(NetworkODE, y0, Params['tvec'],dt,  Params)

# Run the actual sim
Offset = 10
Params['Input'][:, 0] = (50 * ((Params['tvec'] >= 0.1) & (Params['tvec'] <= 0.2)) + 60 * ( (Params['tvec'] >= 0.3) & (Params['tvec'] <= 0.4))  + 70 * ((Params['tvec'] >= 0.5) & (Params['tvec'] <= 0.6))) - Offset  # In A
Params['Input'][:, 1] = (60 * ((Params['tvec'] >= 0.1) & (Params['tvec'] <= 0.2)) + 70 * ( (Params['tvec'] >= 0.3) & (Params['tvec'] <= 0.4)) + 60 * ((Params['tvec'] >= 0.5) & (Params['tvec'] <= 0.6))) - Offset  # In A
Params['Input'][:, 2] = (70 * ((Params['tvec'] >= 0.1) & (Params['tvec'] <= 0.2)) + 50 * ( (Params['tvec'] >= 0.3) & (Params['tvec'] <= 0.4)) + 50 * ((Params['tvec'] >= 0.5) & (Params['tvec'] <= 0.6))) - Offset  # In A
y = euler_solver(NetworkODE,  y[-1, :], Params['tvec'],dt, Params)


# Create a new figure
fig, axarr = plt.subplots(2, 1, sharex=True)  # This creates a 2x1 grid of subplots and shares their x-axes

# Plot on the first axis
axarr[0].plot(Params['tvec'], y[:, [8, 12, 16, 20]])  # Python uses 0-based indexing
axarr[0].set_xlim(tspan)
axarr[0].set_ylabel(r'$V_{\mathrm{mem}} \, (\mathrm{U_T})$')  # LaTeX formatted label
axarr[0].legend([r'$1$', r'$2$', r'$3$', r'$4 \, (\mathrm{Inhibit})$'])  # LaTeX formatted legend
FigFormat(1)

# Plot on the second axis
axarr[1].plot(Params['tvec'], Params['Input'])
axarr[1].set_ylabel(r'$I_{\mathrm{ext}} \, (\mathrm{pA})$')  # LaTeX formatted label
axarr[1].set_xlabel(r'$\mathrm{Time} \, (\mathrm{s})$')  # LaTeX formatted xlabel
FigFormat(1)
axarr[1].autoscale(enable=True, axis='x', tight=True)
axarr[1].set_xlim(tspan)

# Display the plots
plt.tight_layout()
plt.show()