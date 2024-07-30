%cannot have clear or clc in first line, becaus doesn't make it a function,
%makes it a script instead
function [powerCurve] = Nathan_Delos_Santos_power_func(a,x)
    a_0 = a(1);
    n = a(2);

    powerCurve = a_0 * (1 - x/0.038) .^(1/n);
end