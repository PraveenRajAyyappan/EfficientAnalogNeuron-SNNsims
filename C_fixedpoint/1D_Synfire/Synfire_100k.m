%% Initialize
clear, clc;

%% Define some constants
Neuronpop =100e3;
dt = 1e-3;
SimTime = 235;

%% Compute derived constants and preallocate data structures
MaxSize = 10*ceil(SimTime/dt);
Counter = 1;
Timeidx = zeros(MaxSize,1);
SpikeIdx = zeros(MaxSize,1);

%% Load in file into struct
fid = fopen('Synfire_100k.txt');
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
scatter(Timeidx*dt,SpikeIdx,'k.');
xlabel('Time (s)');
ylabel('Neuron');
FigFormat(2);

%% Make raster plot
figure;
%subplot(4,1,1:3), scatter(Timeidx*dt,SpikeIdx,'.');
subplot(4,1,1:3), plot([0.071 226.96;226.871 234.998],[1 1;100000,3583],'LineWidth',2)
ylabel('Neuron');
yticks(0:1e4:1e5)
FigFormat(2);
xlim([0 235]);
timevec = 0:max(Timeidx);
Current = zeros(size(timevec));
Current((timevec>=70) & (timevec<=160)) = 80;
subplot(4,1,4), plot(timevec*dt,Current,'linewidth',2)
xlabel('Time (s)');
ylabel('I_{ext} (pA)')
FigFormat(2);
xlim([0 235]);

%% Zoom 1
Scale = [8 6 4 4];
figure;
scatter(Timeidx*dt,SpikeIdx,'.');
FigFormat(2);
xlim([75 175]*dt)
ylim([0 45])
set(gcf, 'Units','centimeters', 'Position',Scale)

%% Zoom 2
Scale = [8 6 4.5 4.5];
figure;
Temp1 = Timeidx*dt;
idx = (Temp1 >= 226.83) & (Temp1 <= 226.96);
Temp1 = Temp1(idx);
Temp2 = SpikeIdx(idx);
scatter(Temp1,Temp2,'.');
FigFormat(2);
xlim([226.83 226.96])
xticks([226.85,226.95])
ylim([9.994e4 10e4])
set(gcf, 'Units','centimeters', 'Position',Scale)