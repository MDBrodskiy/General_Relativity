e0 = 13.6;          
me = 511000;         
alpha = 1 / 137;    
nb = 10e11;        

alpha_2 = @(x) 9.78 * (alpha^2) / (me^2) .* sqrt(x) .* log(x);

beta = @(x) alpha_2(x) .* ((me * e0) / (2 * pi * x)).^(3/2) .* exp(-x);

dXe_dx = @(x, Xe) (3 / (2 * x.^2)) .* ((1 - Xe) .* beta(x) - Xe.^2 .* nb .* alpha_2(x));

x_start = 1; 
x_end = 1000; 
X0 = 1;

[x, Xe] = ode45(dXe_dx, [x_start:1:x_end], X0);

plot(x, Xe, 'b');
xlabel('x = ϵ0 / T');
ylabel('χe');
title('χe versus x');