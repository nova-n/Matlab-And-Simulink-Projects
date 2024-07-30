clc;
clear;
close all;

%%Problem 1
%transfer function: G(s) = (4s+12)/(s+5)
%numerator (output): (4s+12)
%denominator (input): (s+5)
figure;
G = tf([4,12],[1,5]);
bode(G)
legendStuff = ["$G(s) = \frac{4s+12}{s+5}$"];
lgnd= legend(legendStuff);
set(lgnd, 'Interpreter','latex')
lgnd.Location = 'best';
hold off;

%%Problem 2
%transfer function: G(s) = (1)/(13s^2 +2s + 208)
%numerator (output): (1)
%denominator (input): (13s^2 +2s + 208)
figure;
G = tf([1],[13,2,208]);
bode(G)
legendStuff = ["$G(s) = \frac{1}{13s^2 +2s + 208}$"];
lgnd= legend(legendStuff);
set(lgnd, 'Interpreter','latex')
lgnd.Location = 'best';
hold off;

%%Problem 3
%transfer function: G(s) = (1)/((7.493/3)s+0.327)
%numerator (output): (1)
%denominator (input): ((7.493/3)s+0.327)
figure;
G = tf([1],[7.493/3,0.327]);
bode(G) 
text(1/500,-40,"bandwidth = " + num2str(bandwidth(G)) + "rad/s")
%%might be off by 1.6x, because maybe didn't read graph too accurately
legendStuff = ["$G(s) = \frac{1}{\frac{7.493}{3}s+0.327}$"];
lgnd= legend(legendStuff);
set(lgnd, 'Interpreter','latex')
lgnd.Location = 'best';
hold off;

