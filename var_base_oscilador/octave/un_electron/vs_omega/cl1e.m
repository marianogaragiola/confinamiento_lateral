clear all;
close all;

% algunas constantes utiles
alpha = 658.4092645439; a0 = 0.0529177210; eV = 27.21138564;

% Parametros para el problema
N_cuad = 200;
N_base = 100;

sigma   = 20; % en nm
me      = 0.063; % masa efectiva del electron
B_campo = 80; % campo magnetico en tesla
omega   = 0.5; % parametro variacional no lineal
V0      = 0.05; % Profundidad del potencial en eV

omega_vec = linspace(0., 1, 100);

% Parametros dependientes
b = 0.5*(1/sigma)^2;
l = sqrt(2*alpha/B_campo);

% empiezo con los calculos de las matrices
[x, w] = GaussHermite_2(N_cuad);

file = fopen('./resultados/auval.dat', 'w');
fprintf(file, '# autovalores problema \n');
fclose(file);

for omega = omega_vec

  H = Hamiltoniano(N_base, me, omega, sigma, V0, B_campo, x, w);

  e = eig(H);
  e = sort(e);
  e(:) = eV*e(:);

  e = [omega, e(:)'];

  save('-ascii', '-append', './resultados/auval.dat', 'e')

end
