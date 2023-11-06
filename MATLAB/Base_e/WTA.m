%% This Script Simulates the response of a Farquhar-Hasler Neuron in a WTA
% Written by S. Bhattacharyya, 12/25/22 - Rampur, Howrah District, WB India
%% Initialize
clear, clc;
addpath(genpath(pwd));
load('Params.mat');

%% Define Network Connectivity Via the Synapse Weight Matrix 
% Note: Odd number columns are excitatory, even number columns are inhibitory
Params.NeuronPopulation = 4;    % Number of neurons in the network
Ex = 0.26;                      % Excitatory synapse "strength" (Lower value means greater weightage)
Ih = 0.20;                      % Inhibitory synapse "strength" (Lower value means greater weightage)
Weights = [5  5 5  5 5  5 5 Ih;
           5  5 5  5 5  5 5 Ih;
           5  5 5  5 5  5 5 Ih;
           Ex 5 Ex 5 Ex 5 5 5];              % Controls neuron connectivity
Params.ExternalCurrentIdx = [1,2,3];         % External current will be channeled into neurons with these indices 
Params.Synapse_Mat = Params.M_Syn*Ithp*exp(Kp*(VDD-VT0p)-VDD).*...
    exp(-Kp*Weights/UT).*(+(Weights < 2.5)); % Computes the final transfer matrix

%% Estimate y0 by producing an initial guess
y0 = []; 
for i = 1:Params.NeuronPopulation
    y0 = [y0,[100 100]];
end
for i = 1:Params.NeuronPopulation
    y0 = [y0,Params.DC_OFF];
end

%% Define Simulation Params and Run
dt = 0.0001;     % In s
tspan = [0 0.8]; % In s
Params.tvec  = (tspan(1):dt:tspan(2)).';

% Run a quick transient sim with zero input for the true initial conditions 
Params.Input = zeros(length(Params.tvec),3);                               % In pA
[~,y] = ode15s(@(t,y)NetworkODE(t,y,Params),tspan,y0.');

% Run the actual sim
Offset = 10;
Params.Input(:,1) = (50*(+((Params.tvec >= 0.1) & (Params.tvec <= 0.2))) + 60*(+((Params.tvec >= 0.3) & (Params.tvec <= 0.4))) + 70*(+((Params.tvec >= 0.5) & (Params.tvec <= 0.6)))).'-Offset; % In A
Params.Input(:,2) = (60*(+((Params.tvec >= 0.1) & (Params.tvec <= 0.2))) + 70*(+((Params.tvec >= 0.3) & (Params.tvec <= 0.4))) + 60*(+((Params.tvec >= 0.5) & (Params.tvec <= 0.6)))).'-Offset; % In A
Params.Input(:,3) = (70*(+((Params.tvec >= 0.1) & (Params.tvec <= 0.2))) + 50*(+((Params.tvec >= 0.3) & (Params.tvec <= 0.4))) + 50*(+((Params.tvec >= 0.5) & (Params.tvec <= 0.6)))).'-Offset; % In A
[t,y] = ode15s(@(t,y)NetworkODE(t,y,Params),tspan,y(end,:).');

%% Plot Output
figure, hold on;
ax(1) = subplot(2,1,1);
plot(t,y(:,[9,13,17,21]));
axis tight;
xlim(tspan);
ylabel('V_{mem} (U_T)');
legend('1','2','3','4 (Inhibit)')
ax(2) = subplot(2,1,2);
hold on;
plot(Params.tvec,Params.Input);
ylabel('I_{ext} (pA)');
xlabel('Time (s)');
linkaxes(ax,'x');
axis tight;
xlim(tspan);