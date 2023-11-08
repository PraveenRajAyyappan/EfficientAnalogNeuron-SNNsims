clear, clc;

idx = 30:1:300;

for i = 1:length(idx)
    data = importdata("Vmem_float_Iin_"+num2str(idx(i))+".txt");
    [peak_amplitudes, peak_indices] = findpeaks(data(:,1),'MinPeakProminence',1);
    num_peaks(i)= length(peak_amplitudes);
    firingrate(i) = 1/(median(diff(peak_indices))*0.001e-3);
        
end

%Plot the IF curve
figure, plot(idx,firingrate,'k-o')
xlabel('Input Current (pA)');
ylabel('Spike Frequency (Hz)');

