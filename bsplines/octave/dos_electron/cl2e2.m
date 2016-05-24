
clear all
close all


% algunas constantes utiles
alpha = 658.4092645439; a0 = 0.0529177210; eV = 27.21138564;

%% Parametros fisicos del problema
sigma   = 20; % en nm
me      = 0.063; % masa efectiva del electron
B_campo = 80; % campo magnetico en tesla
V0      = 0.05; % Profundidad del potencial en eV

B_campo_vec = 80; %linspace(10, 100, 90);

% Parametros dependientes
b = 0.5*(1/sigma)^2;

%% Parametros para los B-splines
num_intervalos   = 50; % numero de subintervalos
num_break_points = num_intervalos + 1; % numero de break points contando los extremos
kord             = 5; % orden de los bsplines, el grado de los polinomios es kord-1
regularity       = kord-2; % repeticion de los knots, simplemente degenerados por eso kord-2
xMax             = 150; % extremo superior del intervalo
xMin             = -xMax; % extremo inferior del intervalo
N_cuad           = 100;
N_base = num_intervalos + kord - 1;

name = sprintf('./resultados/2e-E_vs_B-2.dat');

file = fopen(name, 'w');
fprintf(file, '# Num_interavaos = %i \n', num_intervalos)
fprintf(file, '# kord = %i \n', kord)
fprintf(file, '# N_base = %i , N_cuad = %i \n', N_base, N_cuad);
fprintf(file, '# sigma = %f \n', sigma);
fprintf(file, '# me = %f \n', me);
fprintf(file, '# B_campo = %f \n', B_campo);
fprintf(file, '# V0 = %f \n', V0);
fprintf(file, '# Intervalo de intergracion [%f, %f] \n', xMin, xMax)
fprintf(file, '# autovalores del problema \n');
fclose(file);

%% Paso a unidades atomicas todo
sigma = sigma/a0; xMax = xMax/a0; xMin = xMin/a0; b = b*a0^2;
V0 = V0/eV;

%% armo la matriz de control points, no queremos el primer y ultimo bsplines ya que
%% queremos que se cumplan las condiciones de contorno
c = zeros(N_base) + diag([0, ones(1, N_base-2), 0]);

%% armamos los knots
[knots, zeta] = kntuniform (num_break_points, kord-1, regularity);

[x, w] = GaussLegendre_2(N_cuad);
x = x'; w = w';

%% Pasamos las variables al rango de integracion
knots = (xMax - xMin)*knots + xMin; % knots estaba en [0, 1]
zeta = (xMax - xMin)*zeta + xMin;   % zeta estaba en [0, 1]

x = 0.5*(xMax - xMin)*x + 0.5*(xMax + xMin); % x estaba en [-1, 1]

%% calculo los bsplines en el vector x
bs = bspeval(kord-1, c, knots, x);

%% para calcular las derivadas
[dc, dknots] = bspderiv(kord-1, c, knots);

dbs = bspeval(kord-2, dc, dknots, x);

%% Ahora paso a calcular las matrices del problema
T = energia_cinetica(dbs, me, w);

%% Calculamos la energia potencial
[V, S] = energia_potencial(bs, sigma, x, w);
S = S(2:N_base-1,2:N_base-1);

N_dim = (N_base-2)^2;

lambda = 0.45;
V_ef = zeros(N_dim);

for B_campo = B_campo_vec

  l = sqrt(2*alpha/B_campo);
  l = l/a0;

  %% Armo el hamitoniano
  H = zeros(N_base-2);
  H = T(2:N_base-1,2:N_base-1) - V0*V(2:N_base-1,2:N_base-1)/(1 + l^2*b);

  H_sim = kron(H, S) + kron(S, H);
  S_sim = kron(S, S);

  % V_ef = Interaccion(bs, x, w, l);

  % [H_sim, S_sim] = simetrizacion(H, S, V_ef, lambda);

  tic();
  [auvec, e] = eigs(H_sim, S_sim, N_dim);
  [e, perm] = sort(diag(e));
  e = eV*e;
  elapsed_time = toc()
  auval = [e'];

  % auvec = auvec(:, perm);

  save('-ascii', '-append', name, 'auval')

end
