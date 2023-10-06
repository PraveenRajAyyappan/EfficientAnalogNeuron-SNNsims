%% This Script Simulates the response of a Farquhar-Hasler Neuron in Synfire Chain
% Written by S. Bhattacharyya (12/28/22 at Howrah, Howrah District, WB India)
%% Initialize
clear, clc;
addpath(genpath(pwd));
load('Params.mat');

%% Define Network Connectivity Via the Synapse Weight Matrix 
% Note: Odd number columns are excitatory, even number columns are inhibitory
Params.NeuronPopulation = 3;
Syn = 0.29; % Between 0.29 and 0.31 is interesting...
Weights = [2.5 2.5 2.5 2.5 Syn 2.5;
           Syn 2.5 2.5 2.5 2.5 2.5;
           2.5 2.5 Syn 2.5 2.5 2.5]; % Controls neuron connectivity
Params.ExternalCurrentIdx = 1;       % External current will be channeled into neurons with these indices 
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

%% Define Simulation Params and Run
dt = 0.000005;     % In s
tspan = [0 0.5]; % In s
Params.tvec  = (tspan(1):dt:tspan(2)).';

% Run a quick transient sim with zero input for the true initial conditions 
Params.Input = zeros(length(Params.tvec),1);                          % In pA
[~,y] = ode15s(@(t,y)NetworkODE(t,y,Params),tspan,y0.');

% Run the actual sim
Params.Input = (20*(+((Params.tvec >= 0.2) & (Params.tvec <= 0.3)))); % In pA
%[t,y] = ode15s(@(t,y)NetworkODE(t,y,Params),tspan,y(end,:).');
SteadyState = y(end,:);
y = zeros(length(Params.Input),size(y,2));
y(1,:) = SteadyState;
% for i = 2:length(Params.Input) % Implicit Euler
%     y(i,:) = fsolve(@(x)EqSolve(x,y(i-1,:).',dt,Params.tvec(i-1),Params),y(i-1,:).');
% end
for i = 2:length(Params.Input) % Forward Euler
    y(i,:) = y(i-1,:) + dt*NetworkODE(Params.tvec(i-1),y(i-1,:).',Params).';
end
t = Params.tvec;

%% Plot Output
figure, hold on;
ax(1) = subplot(2,1,1);
plot(t,y(:,[7,11,15]));
axis tight;
xlim(tspan);
ylabel('V_{mem} (U_Tln(2))');
legend('1','2','3')
ax(2) = subplot(2,1,2);
hold on;
plot(Params.tvec,Params.Input);
ylabel('I_{ext} (pA)');
xlabel('Time (s)');
linkaxes(ax,'x');
axis tight;
xlim(tspan);