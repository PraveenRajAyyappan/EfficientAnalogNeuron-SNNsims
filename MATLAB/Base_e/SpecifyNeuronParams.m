%% This Script Specifies the Parameters for a Farquhar-Hasler Neuron [Run this if you want to change model parameters!!]
% Written by S. Bhattacharyya (12/25/22 at Rampur, Howrah District, WB India)

%% Initialize
clear, clc;

%% Define Neuron (Physical) Params
% Neuron
M = 8;
Params.CMEM = M*100E-15;
Params.CZ   = M*15E-15; 
Params.CNA  = 8*Params.CZ;
Params.CK   = M*50e-15;
Params.CW   = M*1e-15;
Params.CL   = M*1e-15;
Params.M_Na = 1/0.6; % Microns
Params.M_K  = 1/0.6; % Microns
Params.M_Tn = 1/0.6; % Microns
Params.M_Tm = 12/0.6;% Microns
Params.M_Amp= 12/1;  % Microns
Params.M_Th = 1/0.6; % Microns
Params.M_Tl = 1/0.6; % Microns

% Triangle Gen and Synapse
Params.C_Int = 20e-15;  % Triangle Generator Cap
Params.M_SC = 2/0.6;    % Microns
Params.M_D = 1/0.6;     % Microns
Params.M_Syn = 1/0.6;   % Microns

%% Define Voltage Params
% Misc
UT  = 25.0e-3;                      % Thermal voltage
Params.DC_OFF = [46.03 14 26.558 26.558]; % Estimated initial DC Offsets (in UT). Does not need to be perfect.

% Neuron
Params.ENA  = 52;  
Params.VAMP = 57.6;      
Params.VTH  = -10.8;
Params.VTM  = 28.5;
Params.VSAT = 15.4;
Params.VTN  = -27.2;
Params.VGK  = 14;       
Params.EK   = 46;  
Params.VTL  = 18;      
Params.EL   = 46;

% Triangle and Synapse
Params.V_BP_Ex = 2.18/UT; % P bias, excitatory
Params.V_BN_Ex = 0.13/UT; % N bias, excitatory
Params.V_BP_Ih = 2.225/UT;  % P bias, inhibitory
Params.V_BN_Ih = 0.18/UT; % N bias, inhibitory
Params.SpikeThresh = 49;  % In units of UT

%% Transfer Matrices
% Transistor Parameters
VDD  = 2.5/UT;
Kp   = 0.75;
Ithp = 125e-9;
VT0p = 0.7/UT;
Kn   = 0.67;
Ithn = 300e-9;
VT0n = 0.5/UT;

% Define capacitance elements
C_A_Square = (Params.CNA+Params.CZ+Params.CW)*(Params.CZ+Params.CL) - Params.CZ^2;
C_B_Cube = (Params.CMEM+Params.CNA)*C_A_Square - (Params.CNA^2)*(Params.CZ+Params.CL);
C_11 = C_A_Square/C_B_Cube;
C_14 = Params.CNA*Params.CZ/C_B_Cube;
C_18 = Params.CNA*Params.CL/C_B_Cube;
C_27 = C_A_Square/C_B_Cube + 1/Params.CK;
C_31 = Params.CNA*(Params.CZ+Params.CL)/C_B_Cube;
C_34 = Params.CZ/C_A_Square + Params.CNA^2*Params.CZ*(Params.CZ+Params.CL)/(C_A_Square*C_B_Cube);
C_38 = Params.CL/C_A_Square + Params.CNA^2*Params.CL*(Params.CZ+Params.CL)/(C_A_Square*C_B_Cube);
C_44 = (Params.CNA^2*Params.CZ^2)/(C_A_Square*C_B_Cube) + (Params.CZ^2)/(C_A_Square*(Params.CZ+Params.CL)) + 1/(Params.CZ+Params.CL);
C_48 = (Params.CL*Params.CZ-C_A_Square)/(C_A_Square*(Params.CZ+Params.CL))+(Params.CNA^2*Params.CL*Params.CZ)/(C_A_Square*C_B_Cube);

% Compute Neuron Transfer Matrix
Params.Neuron_Mat = (1/UT)*[C_11,C_11,C_11,C_14,-C_14,-C_11,-C_11,-C_18;
                C_11,C_11,C_11,C_14,-C_14,-C_11,-C_27,-C_18;
                C_31,C_31,C_31,C_34,-C_34,-C_31,-C_31,C_38;
                C_14,C_14,C_14,C_44,-C_44,-C_14,-C_14,-C_48].*...
    repmat([1e-12, Params.M_Tl*Ithp*exp(Kp*(VDD-Params.VTL-VT0p) - VDD + Params.EL),...
    Params.M_Na*Ithp*exp(Kp*(VDD-VT0p)-VDD+Params.ENA), Params.M_Amp*Ithp*exp(Kp*(VDD-VT0p)-VDD+Params.VAMP),...
    Params.M_Tm*Ithn*exp(Kn*(Params.VTM-VT0n) - Params.VSAT), Params.M_K*Ithp*exp(Kp*(VDD-VT0p)-VDD),...
    -Params.M_Tn*Ithp*exp(Kp*(VDD-VT0p-Params.VTN)-VDD+Params.VGK), Params.M_Th*Ithp*exp(Kp*(VDD-Params.VTH-VT0p) - VDD)],[4,1]);
Params.EXPMEL = exp(-Params.EL);
Params.EXPMVGK = exp(-Params.VGK);
Params.EXPVSAT = exp(Params.VSAT);

% Triangle Generator Transfer Matrix
Params.Triangle_Mat = (1/(Params.C_Int*UT))*[-Params.M_D*Ithn*exp(-Kn*(VT0n-Params.V_BN_Ex)),-Params.M_SC*Ithp*exp(Kp*(VDD-VT0p-Params.V_BP_Ex));...
    -Params.M_D*Ithn*exp(-Kn*(VT0n-Params.V_BN_Ih)),-Params.M_SC*Ithp*exp(Kp*(VDD-VT0p-Params.V_BP_Ih))];

%% Save Parameter Matrix
save('Params.mat');