%% This Script Simulates the response of a Single Farquhar-Hasler Neuron and Synapse
% Written by S. Bhattacharyya, 12/25/22 - Rampur, Howrah District, WB India
%% Initialize
clear, clc;
addpath(genpath(pwd));
load('Params.mat');

%% Define Network Connectivity Via the Synapse Weight Matrix 
% Note: Odd number columns are excitatory, even number columns are inhibitory
Params.NeuronPopulation = 1;  % Just one neuron
Weights = [2.5 2.5];          % Controls neuron connectivity. Lower values are more weighted
Params.ExternalCurrentIdx = 1;% External current will be channeled into neurons with these indices 
Params.Synapse_Mat = Params.M_Syn*Ithp*exp(Kp*(VDD-VT0p)-VDD).*...
    exp(-Kp*Weights/UT).*(+(Weights < 2.5)); % Computes the final transfer matrix

%% Estimate y0
y0 = []; 
for i = 1:Params.NeuronPopulation
    y0 = [y0,[100 100]];
end
for i = 1:Params.NeuronPopulation
    y0 = [y0,Params.DC_OFF];
end

%% Define Simulation Params and Run (Leak Channel)
dt = 0.0001;     % In s
tspan = [0 1];   % In s
Params.tvec  = (tspan(1):dt:tspan(2)).';
Params.Input = (40*(+((Params.tvec >= 0.040) & (Params.tvec <= 0.5)))); % In pA

% Find Initial Conditions
options = optimoptions('fsolve');
fun = @(y)NetworkODE(0,y,Params);
x0 = fsolve(fun,y0.',options);

% Run the actual sim
[t,y] = ode15s(@(t,y)NetworkODE(t,y,Params),tspan,x0.');

%% Plot Na Clamp Data
% Plot
figure, hold on;
ax(1) = subplot(3,1,1);
plot(t,y(:,3));
axis tight;
xlim(tspan);
ylabel('V_{mem} (U_T)');

ax(2) = subplot(3,1,2);
plot(t,y(:,1:2));
axis tight;
xlim(tspan);
ylabel('Triangle Out (U_T)');
legend('Excite','Inhibit');

ax(3) = subplot(3,1,3);
hold on;
plot(Params.tvec,Params.Input);
ylabel('I_{ext} (pA)');
xlabel('Time (s)');
linkaxes(ax,'x');
axis tight;
xlim(tspan);

%% Plot Channel Current
GateVolt_E = y(:,1)*0.05 + 90; % Feel free to play with this scaling
GateVolt_I = y(:,2)*0.05 + 10; % Feel free to play with this scaling
Current_E = Params.M_Syn*Ithp*exp(Kp*(VDD-VT0p)).*exp(-Kp*(0.3)/UT).*exp(-Kp*0.8*GateVolt_E);
Current_I = Params.M_Syn*Ithp*exp(Kp*(VDD-VT0p)-VDD).*exp(-Kp*(0.3)/UT).*exp(-Kp*0.8*GateVolt_I).*(1-exp(y(:,3)));
figure, hold on;
plot(t,Current_E*1e12);
plot(t,Current_I*1e12);
ylabel('I_{syn} (pA)');
xlabel('Time (s)');
legend('Excitatory','Inhibitory');