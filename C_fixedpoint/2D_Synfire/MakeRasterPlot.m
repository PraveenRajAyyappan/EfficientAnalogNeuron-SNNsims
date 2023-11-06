%% Initialize
clear, clc;
addpath(genpath(pwd))

%% Define some constants
Neuronpop =100e3;
dt = 1e-3;
SimTime = 5;

%% Compute derived constants and preallocate data structures
MaxSize = 10*ceil(SimTime/dt);
Counter = 1;
Timeidx = zeros(MaxSize,1);
SpikeIdx = zeros(MaxSize,1);

%% Load in file into struct
fid = fopen('Synfire_9.txt');
%fid = fopen('Synfire_10.txt');
%fid = fopen('Synfire_11.txt');
%fid = fopen('Synfire_12.txt');
%fid = fopen('Synfire_13.txt');
%fid = fopen('Synfire_15.txt');
MyTitle = 'n_0=100, \sigma_0=1 ms';
%MyTitle = 'n_0=100, \sigma_0=4.5 ms';
%MyTitle = 'n_0=45, \sigma_0=1 ms';
%MyTitle = 'n_0=100, \sigma_0=4.2 ms';
%MyTitle = 'n_0=48, \sigma_0=1 ms';
%MyTitle = 'n_0=100, \sigma_0=4 ms';
tline = fgetl(fid);
i = 1;
NumDataPoints = 0;
while (ischar(tline))
    tline = fgetl(fid);
    if ((~isempty(tline)) && (i < ceil(SimTime/dt)))
        % Populate the data structures
        temp = str2double(split(tline,' ')); 
        NeuronsNowSpiking = temp(~isnan(temp));
        Increment = length(NeuronsNowSpiking);
        Timeidx((NumDataPoints+1):(NumDataPoints + Increment)) = i*ones(Increment,1);
        SpikeIdx((NumDataPoints+1):(NumDataPoints + Increment)) = NeuronsNowSpiking; 
        NumDataPoints = NumDataPoints + Increment;
    end
    i = i + 1;
end
fclose(fid);
Timeidx = Timeidx(1:NumDataPoints);
SpikeIdx = SpikeIdx(1:NumDataPoints);
[UniqueSpikeIdx,~,~] = unique(SpikeIdx);
for i = 1:length(UniqueSpikeIdx)
    CorrespondingTimeIndices = sort(Timeidx(SpikeIdx == UniqueSpikeIdx(i)));
    TimeLogical = [1 == 1; diff(CorrespondingTimeIndices) > 1];
    CorrespondingTimeIndices = CorrespondingTimeIndices(TimeLogical);

    % Generate mask in original vector and nan them
    Mask = (SpikeIdx == UniqueSpikeIdx(i)) & (~ismember(Timeidx,CorrespondingTimeIndices));
    Timeidx(Mask) = NaN;
    SpikeIdx(Mask) = NaN;
end
% Timeidx = Timeidx(idx);
% tvec = ((0:dt:((Nsteps*dt)-dt)))-0.06;
Timeidx = Timeidx(~isnan(Timeidx));
SpikeIdx = SpikeIdx(~isnan(SpikeIdx));

%% Make raster plot
figure, hold on;
%scatter(1000*Timeidx*dt,SpikeIdx,'k.');
scatter(Timeidx*dt,SpikeIdx,'k.');
%xlabel('Time (ms)');
xlabel('Time (s)');
ylabel('Neuron');
title(MyTitle);
FigFormat(2);
%ylim([0 600])
ylim([0 100000])
%xlim(1000*[0.04, 0.08])
xlim([0,3.5])
xticks(0:0.5:3.5)
yticks(0:1e4:1e5)
%set(gcf, 'Units','centimeters', 'Position',[8 6 5 10])
set(gcf, 'Units','centimeters', 'Position',[8 6 10 10])