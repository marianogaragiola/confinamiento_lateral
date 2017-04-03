clear all
close all

l = 1;
delta = 0.5;
%eta = 0.08;

x = linspace(-20, 20, 201);

I = (abs(x)<= 5);

v_ef = zeros(1, length(x));

v_ef(I) = sqrt(0.5*pi)/l*exp(x(I).*x(I)).*(1 - erf(abs(x(I))));

v_ef(~I) = sqrt(0.5)/l*(1./abs(x(~I)) - 0.5./abs(x(~I)).^3);

v = 1./sqrt(x.^2 + delta^2);

idx = find(x==0);

eta = v_ef(idx)/v(idx);

v = eta*v;

% plot(x, v_ef, 'o-', x, v, 'o-')

salida = [x', v', v_ef'];

save('-ascii', 'l=1.dat', 'salida')

