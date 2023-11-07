% Initialize
clear, clc;

% Load in data
Data = readmatrix("Output.txt");
t = 0:0.001:((length(Data)-1)*0.001);

% Plot
figure, hold on;
plot(t,Data(:,1))
xlabel('Time (ms)')
ylabel('V_{mem} U_Tln(2)')