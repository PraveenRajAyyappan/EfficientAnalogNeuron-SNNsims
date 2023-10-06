function dydt = Neuron(t,y,Params,I_in)
% FH Type HH Neuron ODE
% Written by S. Bhattacharyya (12/28/22 at Howrah, Howrah District, WB India)
% y order: v_mem, v_k, v_g, v_na
Kp    = 0.75;
I_Vec = [I_in;
    1 - Params.EXPMEL*MyPwrFunc(y(1));
    MyPwrFunc(-Kp*y(4));                 % This is an approximation (dropped a -exp(y(1)-Kp*y(4))term)
    MyPwrFunc(-Kp*y(3));                 % This is an approximation (dropped a -exp(y(4)-Kp*y(3)) term)
    1 - Params.EXPVSAT*MyPwrFunc(-y(4));
    MyPwrFunc(y(1)-Kp*(y(2)));           % This is an approximation (dropped a -exp(Params.EK-Kp*y(2)) term)
    1 - Params.EXPMVGK*MyPwrFunc(y(2));
    MyPwrFunc(y(3))-MyPwrFunc(y(4))];
dydt = Params.Neuron_Mat*I_Vec;
end