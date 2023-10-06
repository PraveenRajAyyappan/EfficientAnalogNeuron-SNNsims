function dydt = Triangle(t,y,Params,Vmem)
% Triangle Generator ODE (Randall Alpha)
% Parameters
VDD = 144.2695; % In UTln(2) of course...

% y order: Excitatory, Inhibitory
zeta = +(Vmem >= Params.SpikeThresh);      % Perform Comparison (level crossing detection)
dydt = Params.Triangle_Mat*[zeta;zeta-1];  % Compute actual derivative

% Error checking (Make sure no change that will push the output beyond the rails)
dydt((y >= VDD).*sign(dydt) == 1) = 0;
dydt(-(y <= 0).*sign(dydt) == 1) = 0;
end