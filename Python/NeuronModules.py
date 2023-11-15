# This script contains the functions and parameters necessary to properly run the WTA, Synfire, and OneNeuron scripts.
import numpy as np
import matplotlib.pyplot as plt

def euler_solver(ode_function, y0, t, dt, Params):
    # Initialize the array to hold the solution
    y = np.zeros((len(t), len(y0)))
    y[0] = y0

    # Iterate over all time points
    for i in range(1, len(t)):
        # Calculate the derivatives at the current time point
        dydt = ode_function(t[i - 1], y[i - 1], Params)
        # Update the state variables using the Euler method
        y[i] = y[i - 1] + dydt * dt

    return y

def Triangle(t, y, Params, Vmem):
    # Triangle Generator ODE (Randall Alpha)
    # Parameters
    VDD = 100  # In UT

    # y order: Excitatory, Inhibitory
    zeta = (Vmem >= Params['SpikeThresh'])  # Perform Comparison (level crossing detection)
    dydt = np.dot(Params['Triangle_Mat'], np.array([[zeta], [zeta - 1]]))  # Compute actual derivative

    # Error checking (Make sure no change that will push the output beyond the rails)s
    dydt[((y >= VDD) & (np.sign(dydt) == 1))] = 0
    dydt[((y <= 0) & (np.sign(dydt) == -1))] = 0
    return dydt

def Neuron(t, y, Params, I_in):
    # FH Type HH Neuron ODE
    # y order: v_mem, v_k, v_g, v_na
    Kp = 0.75

    I_Vec = np.array([
        I_in,
        1 - Params['EXPMEL'] * np.exp(y[0,0]),
        np.exp(-Kp * y[3,0]) * (1 - np.exp(y[0,0] - Params['ENA'])),
        np.exp(-Kp * y[2,0]),
        1 - Params['EXPVSAT'] * np.exp(-y[3,0]),
        np.exp(y[0,0] - Kp * y[1,0]) - np.exp(Params['EK'] - Kp * y[1,0]),
        1 - Params['EXPMVGK'] * np.exp(y[1,0]),
        np.exp(y[2,0]) - np.exp(y[3,0])
    ])
    
    dydt = np.dot(Params['Neuron_Mat'], I_Vec).reshape(-1, 1)
    return dydt

def NeuronTest(y, t, Params):
    y = y.reshape(-1, 1)
    InputCurrent = np.interp(t, Params['tvec'], Params['Input'])
    dydt = Neuron(t, y, Params, InputCurrent)
    return dydt.ravel()

def FigFormat(type):
    # Function to format figures into standardized sizes
    # Get current figure handle
    # Note: type==1 --> Portrait mode, type==2 --> Landscape mode
    fig = plt.gcf()
    
    # Set default line widths
    plt.box(True)
    plt.rc('lines', linewidth=1)
    plt.rc('axes', linewidth=1)
    
    # Resize figure (assuming you want to set size in centimeters, converting to inches)
    if type == 1:
        fig.set_size_inches(16/2.54, 16/2.54)  # Converting cm to inches
    else:
        fig.set_size_inches(16.18/2.54, 10/2.54)  # Converting cm to inches
    
    # Set background color to white
    fig.set_facecolor('w')
    
    # Set font sizes
    plt.rc('font', size=12)
    plt.title(plt.gca().get_title(), fontsize=12)

def NetworkODE(t,y, Params):
    y = y.reshape(-1, 1)   # Change the shape of the array to switch the number of columns and rows
    VDD = 100   # What the pfet well is tied to, in the UT 
    Kp = 0.75   # Kappa
    CG_Coupling = 0.8   # Control of the Gate Coupling Factor
    EffCoupling = CG_Coupling * Kp  # Combination of both the CG-gate and gate-channel coupling
    Tri_Top_Off = 90    # Resting offset of the control gate of the excitatory synapse
    Tri_Bot_Off = 10    # Resting offset of the control gate of the inhibitory synapse
    TriScaling = 0.05   # Scaling of the triangle generator outputs

    # Extract and formatting information
    dydt = np.zeros_like(y)     # Preallocation
    Vmem = y[(2*Params['NeuronPopulation'])::4,:]   # The membrane voltage of each neuron

    # Compute input current
    # Finding the indices of excitatory and inhibitory synapses
    idx = np.arange(2*Params['NeuronPopulation']).reshape(-1, 1)
    Bot_Idx = (idx % 2).astype("double")    # Inhibitory
    Top_Idx = 1 - Bot_Idx   # Excitatory
    Triangle_Scaled = TriScaling *y[0:(2*Params['NeuronPopulation']),:] + Tri_Top_Off*Top_Idx + Tri_Bot_Off*Bot_Idx # Scaling of Triangle generator outputs
    ScaleMat = np.zeros_like(Params['Synapse_Mat'])     # Preallocation

    top_idx_cols = idx[Top_Idx == 1]    # Assigns the indices [0 2 4 6] for excitatory synapses
    bot_idx_cols = idx[Bot_Idx == 1]    # Assigns the indices [1 3 5 7] for inhibitory synapses

    # Assign values using integer indexing
    ScaleMat[:, top_idx_cols] = np.exp(VDD) * np.ones((ScaleMat.shape[0], top_idx_cols.size))   # Appending additional factors for excitatory synapses
    ScaleMat[:, bot_idx_cols] = (1 - np.exp(Vmem)).reshape(-1, 1) * np.ones((ScaleMat.shape[0], bot_idx_cols.size))     # Appending additional factors for inhibitory synapses

    # Neuron input current calculation
    InputCurrent = 1e12 * np.dot(np.multiply(Params['Synapse_Mat'] , ScaleMat) , np.exp(-EffCoupling * Triangle_Scaled))    # From neuron-neuron connections
    
    for i, idx in enumerate(Params['ExternalCurrentIdx']):  # Append the contributions of external disturbances
        InputCurrent[idx-1,0] += np.interp(t, Params['tvec'], Params['Input'][:, i])
    
    # Solve for the triangle generator outputs
    for i in range(Params['NeuronPopulation']):
        start_idx = 2*i
        dydt[start_idx:(start_idx+2),:] = Triangle(t, y[start_idx:(start_idx+2),:], Params, Vmem[i,0])  # Solve for dydt of the triangle generators for the current neuron
    
    # Solve for the neuron update
    for i in range(Params['NeuronPopulation']):
        StartIdx = (2 * Params['NeuronPopulation']) + 4 *i  # Corresponds to the index of the first state variable of the current neuron
        dydt[StartIdx:(StartIdx+4),:] = Neuron(t, y[StartIdx:(StartIdx+4),:], Params, InputCurrent[i,0])    # Solve for dydt of the current neuron
    return dydt.ravel()

# Initialize
np.random.seed(0)

# Define Neuron (Physical) Params
Params = {}

# Neuron
M = 8
Params['CMEM'] = M * 100E-15
Params['CZ'] = M * 15E-15 
Params['CNA'] = 8 * Params['CZ']
Params['CK'] = M * 50e-15
Params['CW'] = M * 1e-15
Params['CL'] = M * 1e-15
Params['M_Na'] = 1 / 0.6  # Microns
Params['M_K'] = 1 / 0.6  # Microns
Params['M_Tn'] = 1 / 0.6  # Microns
Params['M_Tm'] = 12 / 0.6  # Microns
Params['M_Amp'] = 12 / 1  # Microns
Params['M_Th'] = 1 / 0.6  # Microns
Params['M_Tl'] = 1 / 0.6  # Microns

# Triangle Gen and Synapse
Params['C_Int'] = 20e-15  # Triangle Generator Cap
Params['M_SC'] = 2 / 0.6  # Microns
Params['M_D'] = 1 / 0.6  # Microns
Params['M_Syn'] = 1 / 0.6  # Microns

# Define Voltage Params
# Misc
UT = 25.0e-3  # Thermal voltage
Params['DC_OFF'] = [46.03, 14, 26.558, 26.558]  # Estimated initial DC Offsets (in UT). Does not need to be perfect.

# Neuron
Params['ENA'] = 52
Params['VAMP'] = 57.6
Params['VTH'] = -10.8
Params['VTM'] = 28.5
Params['VSAT'] = 15.4
Params['VTN'] = -27.2
Params['VGK'] = 14
Params['EK'] = 46
Params['VTL'] = 18
Params['EL'] = 46

# Triangle and Synapse
Params['V_BP_Ex'] = 2.18 / UT  # P bias, excitatory
Params['V_BN_Ex'] = 0.13 / UT  # N bias, excitatory
Params['V_BP_Ih'] = 2.225 / UT  # P bias, inhibitory
Params['V_BN_Ih'] = 0.18 / UT  # N bias, inhibitory
Params['SpikeThresh'] = 49  # In units of UT

# Transfer Matrices
# Transistor Parameters
VDD = 2.5 / UT
Kp = 0.75
Ithp = 125e-9
VT0p = 0.7 / UT
Kn = 0.67
Ithn = 300e-9
VT0n = 0.5 / UT

# Define capacitance elements
C_A_Square = (Params['CNA'] + Params['CZ'] + Params['CW']) * (Params['CZ'] + Params['CL']) - Params['CZ']**2
C_B_Cube = (Params['CMEM'] + Params['CNA']) * C_A_Square - (Params['CNA']**2) * (Params['CZ'] + Params['CL'])
C_11 = C_A_Square / C_B_Cube
C_14 = Params['CNA'] * Params['CZ'] / C_B_Cube
C_18 = Params['CNA'] * Params['CL'] / C_B_Cube
C_27 = C_A_Square / C_B_Cube + 1 / Params['CK']
C_31 = Params['CNA'] * (Params['CZ'] + Params['CL']) / C_B_Cube
C_34 = Params['CZ'] / C_A_Square + Params['CNA']**2 * Params['CZ'] * (Params['CZ'] + Params['CL']) / (C_A_Square * C_B_Cube)
C_38 = Params['CL'] / C_A_Square + Params['CNA']**2 * Params['CL'] * (Params['CZ'] + Params['CL']) / (C_A_Square * C_B_Cube)
C_44 = (Params['CNA']**2 * Params['CZ']**2) / (C_A_Square * C_B_Cube) + (Params['CZ']**2) / (C_A_Square * (Params['CZ'] + Params['CL'])) + 1 / (Params['CZ'] + Params['CL'])
C_48 = (Params['CL'] * Params['CZ'] - C_A_Square) / (C_A_Square * (Params['CZ'] + Params['CL'])) + (Params['CNA']**2 * Params['CL'] * Params['CZ']) / (C_A_Square * C_B_Cube)

# Compute Neuron Transfer Matrix
Params['Neuron_Mat'] = (1 / UT) * np.array([[C_11, C_11, C_11, C_14, -C_14, -C_11, -C_11, -C_18],
                                             [C_11, C_11, C_11, C_14, -C_14, -C_11, -C_27, -C_18],
                                             [C_31, C_31, C_31, C_34, -C_34, -C_31, -C_31, C_38],
                                             [C_14, C_14, C_14, C_44, -C_44, -C_14, -C_14, -C_48]]) * \
    np.array([1e-12, Params['M_Tl'] * Ithp * np.exp(Kp * (VDD - Params['VTL'] - VT0p) - VDD + Params['EL']),
              Params['M_Na'] * Ithp * np.exp(Kp * (VDD - VT0p) - VDD + Params['ENA']),
              Params['M_Amp'] * Ithp * np.exp(Kp * (VDD - VT0p) - VDD + Params['VAMP']),
              Params['M_Tm'] * Ithn * np.exp(Kn * (Params['VTM'] - VT0n) - Params['VSAT']),
              Params['M_K'] * Ithp * np.exp(Kp * (VDD - VT0p) - VDD),
              -Params['M_Tn'] * Ithp * np.exp(Kp * (VDD - VT0p - Params['VTN']) - VDD + Params['VGK']),
              Params['M_Th'] * Ithp * np.exp(Kp * (VDD - Params['VTH'] - VT0p) - VDD)])

Params['EXPMEL'] = np.exp(-Params['EL'])
Params['EXPMVGK'] = np.exp(-Params['VGK'])
Params['EXPVSAT'] = np.exp(Params['VSAT'])

# Triangle Generator Transfer Matrix
Params['Triangle_Mat'] = (1 / (Params['C_Int'] * UT)) * np.array([[-Params['M_D'] * Ithn * np.exp(-Kn * (VT0n - Params['V_BN_Ex'])),
                                                                  -Params['M_SC'] * Ithp * np.exp(Kp * (VDD - VT0p - Params['V_BP_Ex']))],
                                                                 [-Params['M_D'] * Ithn * np.exp(-Kn * (VT0n - Params['V_BN_Ih'])),
                                                                  -Params['M_SC'] * Ithp * np.exp(Kp * (VDD - VT0p - Params['V_BP_Ih']))]])
