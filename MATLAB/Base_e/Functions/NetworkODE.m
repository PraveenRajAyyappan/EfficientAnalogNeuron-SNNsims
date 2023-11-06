function dydt = NetworkODE(t,y,Params)
% y order: Triangle Inputs, Neuron State
%% Some More Parameters...
VDD         = 100; % In UT of course...
Kp          = 0.75;% Kappa
CG_Coupling = 0.8; % Control Gate Coupling Factor
EffCoupling = CG_Coupling*Kp; % Combination of both the CG-gate and gate-channel coupling
Tri_Top_Off = 90;  % Resting offset of the control gate of the excitatory synapse
Tri_Bot_Off = 10;  % Resting offset of the control gate of the inhibitory synapse
TriScaling  = 0.05;% Scaling of the triangle generator outputs

%% Extract and format stuff
dydt = zeros(length(y),1);                     % Preallocation
Vmem = y((2*Params.NeuronPopulation+1):4:end); % The membrane voltage of each neuron

%% Compute input current
% Finding indices of excitatory and inhibitory synapses
idx = 1:2*Params.NeuronPopulation;
Top_Idx = mod(idx,2); % Excitatory
Bot_Idx = 1-Top_Idx;  % Inhibitory
Triangle_Scaled = y(1:(2*Params.NeuronPopulation))*TriScaling +...
    Top_Idx.'*Tri_Top_Off + Bot_Idx.'*Tri_Bot_Off;                                    % Scaling of Triangle generator outputs
ScaleMat = zeros(size(Params.Synapse_Mat));                                           % Preallocation
ScaleMat(:,idx(Top_Idx == 1)) = repmat(exp(VDD)*ones(size(Vmem)),[1,Params.NeuronPopulation]);% Appending additional factors for excitatory synapses
ScaleMat(:,idx(Bot_Idx == 1)) = repmat(1-exp(Vmem),[1,Params.NeuronPopulation]);              % Appending additional factors for inhibitory synapses

% Neuron input current calculation
InputCurrent = 1e12*(Params.Synapse_Mat.*ScaleMat)*exp(-EffCoupling*Triangle_Scaled); % Due to neuron-neuron connections
for i = 1:length(Params.ExternalCurrentIdx)                                           % Append the contributions of external disturbances
    InputCurrent(Params.ExternalCurrentIdx(i)) = InputCurrent(Params.ExternalCurrentIdx(i)) +...
        interp1(Params.tvec,Params.Input(:,i),t);
end

%% Solve for triangle generator outputs
for i = 1:Params.NeuronPopulation
    try
    dydt((2*i-1):(2*i)) = Triangle(t,y((2*i-1):(2*i)),Params,Vmem(i));     % Solve for dydt of the triangle generators for the current neuron
    catch

    end
end

%% Solve for the neuron update
for i = 1:Params.NeuronPopulation
    StartIdx = (2*Params.NeuronPopulation+1)+4*(i-1); % Corresponds to the index of the first state variable of the current neuron
    dydt(StartIdx:(StartIdx+3)) = Neuron(t,y(StartIdx:(StartIdx+3)),Params,InputCurrent(i)); % Solve for dydt of the current neuron
end
end