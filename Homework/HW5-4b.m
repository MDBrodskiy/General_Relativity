clear all;
close all;

e0 = 13.6;          
me = 511000;         
alpha = 1 / 137;    
nb = 1e-6;        

alpha_2 = @(x) 9.78 * (alpha^2) / (me^2) * sqrt(x) * log(x);

beta = @(x) alpha_2(x) * ((me * e0) / (2 * pi * x))^(3/2) * exp(-x);

%dXe_dx = @(x, Xe) (3 / (2 * x^2)) * ((1 - Xe) * beta(x) - Xe^2 * nb * alpha_2(x));
dXe_dx = @(x, Xe) (1 - Xe) .* (1.12e20) .* sqrt(x) .* log(x) .* exp(-x) - Xe^2 .* (3.99e4) .* x^(-2) .* log(x);

x_start = 1; 
x_end = 10; 
X0 = 1;

[x, Xe] = ode45(dXe_dx, [x_start:x_end], X0);

plot(x, Xe, 'b');
xlabel('x = ϵ0 / T');
ylabel('χe');
title('χe versus a');