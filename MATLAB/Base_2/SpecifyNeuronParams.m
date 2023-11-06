%% This Script Specifies the Parameters for a Farquhar-Hasler Neuron [Run this if you want to change model parameters!!]
% Written by S. Bhattacharyya (12/28/22 at Howrah, Howrah District, WB India)

%% Initialize
clear, clc;

%% Define Neuron (Physical) Params
% Neuron
M = 8;
Params.CMEM = M*100E-15;
Params.CZ   = M*15E-15; 
Params.CNA  = 8*Params.CZ;
Params.CK   = M*50e-15;
Params.M_Na = 1/0.6; % Microns
Params.M_K  = 1/0.6; % Microns
Params.M_Tn = 1/0.6; % Microns
Params.M_Tm = 12/0.6;% Microns
Params.M_Amp= 12/1;  % Microns
Params.M_Th = 1/0.6; % Microns
Params.M_Tl = 1/0.6; % Microns

% Triangle Gen and Synapse
Params.C_Int = 20e-15;  % Triangle Generator Cap
Params.M_SC  = 2/0.6;   % Microns
Params.M_D   = 1/0.6;   % Microns
Params.M_Syn = 1/0.6;   % Microns

%% Define Voltage Params
% Misc
UT  = 25.0e-3;                   % Thermal voltage
ScaleFactor = 1/(log(2)*UT);
OffOff = [15 0 0 0];
Params.DC_OFF = [56.0  19.6  47.6  47.6]+OffOff; % Esimated initial DC Offsets in UTln(2) 
% MAKE SURE element 3 and 4 are the SAME ^^

% Neuron
Params.ENA  = 1.385*ScaleFactor;  
Params.VAMP = 1.535*ScaleFactor;      
Params.VTH  = -0.07*ScaleFactor;
Params.VTM  = 0.740*ScaleFactor;
Params.VSAT = 0.430*ScaleFactor;
Params.VTN  = -0.693*ScaleFactor;
Params.VGK  = 0.341*ScaleFactor;       
Params.EK   = 0.830*ScaleFactor;   
Params.VTL  = 0.340*ScaleFactor;      
Params.EL   = 0.97*ScaleFactor;

% Triangle and Synapse
Params.V_BP_Ex = 2.18*ScaleFactor; % P bias, excitatory
Params.V_BN_Ex = 0.13*ScaleFactor; % N bias, excitatory
Params.V_BP_Ih = 2.225*ScaleFactor;% P bias, inhibitory
Params.V_BN_Ih = 0.18*ScaleFactor; % N bias, inhibitory
Params.SpikeThresh = -3.2;         % In units of UTln(2)

%% Transfer Matrices
% Transistor Parameters
VDD  = 2.5*ScaleFactor;
Kp   = 0.75;
Ithp = 125e-9;
VT0p = 0.7*ScaleFactor;
Kn   = 0.67;
Ithn = 300e-9;
VT0n = 0.5*ScaleFactor;

% Compute Neuron Transfer Matrix
Params.Neuron_Mat = (1/(Params.CMEM*UT*log(2)))*[1, 1, 1, 1, -1, -1, -1, 0;...
    1, 1, 1, 1, -1, -1, -(1+Params.CMEM/Params.CK), 0;...
    1, 1, 1, 1+Params.CMEM/Params.CNA, -(1+Params.CMEM/Params.CNA), -1, -1, 0;...
    1, 1, 1, 1+Params.CMEM/Params.CNA+Params.CMEM/Params.CZ,-(1+Params.CMEM/Params.CNA+Params.CMEM/Params.CZ), -1, -1, Params.CMEM/Params.CZ].*...
    repmat([1e-12, Params.M_Tl*Ithp*2^(Kp*(VDD-Params.VTL-VT0p) - VDD + Params.EL),...
    Params.M_Na*Ithp*2^(Kp*(VDD-VT0p-Params.DC_OFF(4))-VDD+Params.ENA), Params.M_Amp*Ithp*2^(Kp*(VDD-VT0p-Params.DC_OFF(3))-VDD+Params.VAMP),...
    Params.M_Tm*Ithn*2^(Kn*(Params.VTM-VT0n) - Params.VSAT), Params.M_K*Ithp*2^(Kp*(VDD-VT0p-Params.DC_OFF(2))-VDD+Params.DC_OFF(1)),...
    -Params.M_Tn*Ithp*2^(Kp*(VDD-VT0p-Params.VTN)-VDD+Params.VGK), Params.M_Th*Ithp*2^(Kp*(VDD-Params.VTH-VT0p) - VDD + Params.DC_OFF(3))],[4,1]);
Params.EXPMEL = 2^(Params.DC_OFF(1)-Params.EL);
Params.EXPMVGK = 2^(Params.DC_OFF(2)-Params.VGK);
Params.EXPVSAT = 2^(-Params.DC_OFF(4)+Params.VSAT);

% Triangle Generator Transfer Matrix
Params.Triangle_Mat = (1/(Params.C_Int*UT*log(2)))*[-Params.M_D*Ithn*2^(-Kn*(VT0n-Params.V_BN_Ex)),-Params.M_SC*Ithp*2^(Kp*(VDD-VT0p-Params.V_BP_Ex));...
    -Params.M_D*Ithn*2^(-Kn*(VT0n-Params.V_BN_Ih)),-Params.M_SC*Ithp*2^(Kp*(VDD-VT0p-Params.V_BP_Ih))];
Params.Bot_Scale = -2^Params.DC_OFF(1);

%% Save Parameter Matrix
save('Params.mat');