%% This Script Simulates the response of a Single Farquhar-Hasler Neuron and Synapse
% Written by S. Bhattacharyya (12/28/22 at Howrah, Howrah District, WB India)
%% Initialize
clear, clc;
addpath(genpath(pwd));
load('Params.mat');

%% Define Network Connectivity Via the Synapse Weight Matrix 
% Note: Odd number columns are excitatory, even number columns are inhibitory
Params.NeuronPopulation = 1;  % Just one neuron
Weights = [2.5 2.5];          % Controls neuron connectivity. Lower values are more weighted
Params.ExternalCurrentIdx = 1;% External current will be channeled into neurons with these indices 
Params.Synapse_Mat = Params.M_Syn*Ithp*2^(Kp*(VDD-VT0p)-VDD).*...
    2.^(-Kp*Weights/(UT*log(2))).*(+(Weights < 2.5)); % Computes the final transfer matrix

%% Estimate y0
y0 = []; 
for i = 1:Params.NeuronPopulation
    y0 = [y0,[VDD VDD]];
end
for i = 1:Params.NeuronPopulation
    y0 = [y0,-OffOff];
end

%% Define Simulation Params and Run (Leak Channel)
dt = 0.000005;   % In s
tspan = [0 0.2]; % In s
Params.tvec  = (tspan(1):dt:tspan(2)).';

% Run a quick transient sim with zero input for the true initial conditions 
Params.Input = zeros(length(Params.tvec),1);                            % In pA
[~,y] = ode15s(@(t,y)NetworkODE(t,y,Params),tspan,y0.');

% Run the actual sim
Params.Input = (45*(+((Params.tvec >= 0.010) & (Params.tvec <= 0.1)))); % In pA
[t,y] = ode15s(@(t,y)NetworkODE(t,y,Params),tspan,y(end,:).');

%% Plot Data
% Plot
figure, hold on;
ax(1) = subplot(3,1,1);
plot(t,y(:,3));
axis tight;
xlim(tspan);
ylabel('V_{mem} (U_Tln(2), Offset)');

ax(2) = subplot(3,1,2);
plot(t,y(:,1:2));
axis tight;
xlim(tspan);
ylabel('Triangle Out (U_Tln(2))');
legend('Excite','Inhibit');

ax(3) = subplot(3,1,3);
hold on;
plot(Params.tvec,Params.Input);
ylabel('I_{ext} (pA)');
xlabel('Time (s)');
linkaxes(ax,'x');
axis tight;
xlim(tspan);