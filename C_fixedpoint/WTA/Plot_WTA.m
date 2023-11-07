% Initialize
clear, clc;

% Load in data
Data = readmatrix("WTA.txt");
t = 0:0.001:((length(Data)-1)*0.001);

% Plot
figure, hold on;
for i = 1:5
    plot(t,Data(:,i))
end
legend('1','2','3','4','Inhibitory')
xlabel('Time (ms)')
ylabel('V_{mem} U_Tln(2)')