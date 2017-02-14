clear all
close all

l = 10;

x = linspace(-20, 20, 201);

I = (abs(x)<= 5);

v = zeros(1, length(x));

v(I) = sqrt(0.5*pi)/l*exp(x(I).*x(I)).*(1 - erf(abs(x(I))));

v(~I) = sqrt(0.5)/l*(1./abs(x(~I)) - 0.5./abs(x(~I)).^3);

plot(x, v, 'o-')
