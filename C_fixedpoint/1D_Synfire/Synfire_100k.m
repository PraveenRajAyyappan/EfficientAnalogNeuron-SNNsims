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
subplot(4,1,1:3), scatter(Timeidx*dt,SpikeIdx,'.');
ylabel('Neuron');
yticks(0:1e4:1e5)
xlim([0 235]);
timevec = 0:max(Timeidx);
Current = zeros(size(timevec));
Current((timevec>=70) & (timevec<=160)) = 80;
subplot(4,1,4), plot(timevec*dt,Current,'linewidth',2)
xlabel('Time (s)');
ylabel('I_{ext} (pA)')
xlim([0 235]);