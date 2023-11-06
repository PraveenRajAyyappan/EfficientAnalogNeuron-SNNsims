function MyPower = MyPwrFunc(x)
Cfac = 0.25;       % Divide by 4...
Delta = x-floor(x);
MyPower = 2.^floor(x).*(1+Delta).*(Cfac*Delta.^2-Cfac.*Delta+1);
end