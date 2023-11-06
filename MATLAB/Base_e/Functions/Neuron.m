function dydt = Neuron(t,y,Params,I_in)
% FH Type HH Neuron ODE
% y order: v_mem, v_k, v_g, v_na
Kp    = 0.75;
I_Vec = [I_in;
    1 - Params.EXPMEL*exp(y(1));
    exp(-Kp*y(4))*(1-exp(y(1)-Params.ENA));
    exp(-Kp*y(3));            
    1 - Params.EXPVSAT*exp(-y(4));
    exp(y(1)-Kp*(y(2))) - exp(Params.EK-Kp*(y(2)));
    1 - Params.EXPMVGK*exp(y(2));
    exp(y(3))-exp(y(4))];
dydt = Params.Neuron_Mat*I_Vec;
end