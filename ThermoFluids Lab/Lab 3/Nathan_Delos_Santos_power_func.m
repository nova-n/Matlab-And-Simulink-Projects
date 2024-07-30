%cannot have clear or clc in first line, becaus doesn't make it a function,
%makes it a script instead
function [powerCurve] = Nathan_Delos_Santos_power_func(a,x)
    a_0 = a(1);
    a_1 = a(2);
    n = a(3);

    powerCurve = a_0 * a_1*(x).^(1/n);
end