clear all
close all

a0 = 0.0529177210; eV = 27.21138564; c_light = 137.035999074492; ua_to_T = 1.72e3;

me = 0.041;
Rmax = 50.0;
Bi = 10.0;
Bf = 20.0;

bcampo_vec = linspace(Bi, Bf, 6);

r = linspace(0, Rmax, 50);
r = r';

for i = 1:6
  bcampo = bcampo_vec(i);

  name = sprintf('./resultados/landau_level-B_%3.1f.dat', bcampo);
  file = fopen(name, 'w');

  omega = 0.5*bcampo/(me*c_light*ua_to_T);

  lcampo = sqrt(0.5/(me*omega))*a0;

  phi = 1/lcampo^2*exp(-r.^2/(4*lcampo^2));
  phi = r.*phi.^2;

  ll = [r, phi];
  save('-ascii', name, 'll');
end
