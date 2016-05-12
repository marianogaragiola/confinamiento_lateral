close all
clear all

format long

h = 1.0545716E-27; % erg*s
qe = 4.80320427E-10; % Fr
me = 9.10938215E-28; % g
me_red = 0.063; % adimencional
c = 2.99792458E10; % cm/s

erg_to_eV = 6.241506363e+11; % 1 erg <----> 6.241506363e+11 eV
G_to_T = 1E-4; % 1G <----> 1E-4 T
cm_to_nm = 1E7;% 1cm <----> 10000000nm

alpha1 = 2.0*h*c/qe; % G*cm^2

alpha2 = alpha1*G_to_T*cm_to_nm^2; % T*nm^2

beta1 = (h*qe)/(me_red*me*c); % erg/G

beta2 = beta1*erg_to_eV; % eV/G

beta3 = beta2/G_to_T; % eV/T

B = linspace(0, 25, 25)';

salida = [B];

for n = 0:10
  salida = [salida, 1e3*beta3*(n+1/2)*B];
end

save('-ascii', 'landau_levels.dat', 'salida')

% plot(B, 1e3*beta3*B);
