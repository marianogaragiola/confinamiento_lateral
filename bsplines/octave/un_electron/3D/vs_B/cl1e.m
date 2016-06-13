
clear all
close all
tic();

% algunas constantes utiles
alpha = 658.4092645439; a0 = 0.0529177210; eV = 27.21138564;

%% Parametros fisicos del problema
sigma   = 20; % en nm
me      = 0.063; % masa efectiva del electron
B_campo = 80; % campo magnetico en tesla
V0      = 0.05; % Profundidad del potencial en eV

V0_vec = linspace(0, 0.05, 10);

% Parametros dependientes
b = 0.5*(1/sigma)^2;

%% Parametros para los B-splines
num_intervalos   = 10; % numero de subintervalos
num_break_points = num_intervalos + 1; % numero de break points contando los extremos
kord             = 5; % orden de los bsplines, el grado de los polinomios es kord-1
regularity       = kord-2; % repeticion de los knots, simplemente degenerados por eso kord-2
xMax             = 150; % extremo superior del intervalo
xMin             = 0; % extremo inferior del intervalo
N_cuad           = 100;
N_base = num_intervalos + kord - 1;


%% Algunas ctes utiles
h = 1.0545716E-27; % erg*s
qe = 4.80320427E-10; % Fr
melectron = 9.10938215E-28; % g
me_red = 0.063; % adimencional
c = 2.99792458E10; % cm/s
erg_to_eV = 6.241506363e+11; % 1 erg <----> 6.241506363e+11 eV
G_to_T = 1E-4; % 1G <----> 1E-4 T
cm_to_nm = 1E7;% 1cm <----> 10000000nm

%% Frecuencias en cada coordenada
%% Las frecuencias de oscilacion debido al campo magnetico son
%% \omega = e^2 B^2/(4*m*c^2);
omega = qe^2*B_campo.^2/(4*melectron*c*G_to_T.^2); % deberia de estar en ergios/cm^2
omega = omega*erg_to_eV/cm_to_nm^2; % unidades atomicas

omega_x = 0.0;
omega_y = 0.0;
omega_z = 0.0;

name = sprintf('./resultados/1e-E_vs_V0.dat');

file = fopen(name, 'w');
fprintf(file, '# Num_interavaos = %i \n', num_intervalos)
fprintf(file, '# kord = %i \n', kord)
fprintf(file, '# N_base = %i , N_cuad = %i \n', N_base, N_cuad);
fprintf(file, '# sigma = %f \n', sigma);
fprintf(file, '# me = %f \n', me);
fprintf(file, '# Intervalo de intergracion [%f, %f] \n', xMin, xMax)
fprintf(file, '# autovalores del problema \n');
fclose(file);

%% Paso a unidades atomicas todo
sigma = sigma/a0; xMax = xMax/a0; xMin = xMin/a0; b = b*a0^2;
l = sqrt(2*alpha/B_campo); l = l/a0;
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
Tx = energia_cinetica(dbs, me, w);
Tx = Tx(2:N_base-1,2:N_base-1);

%% Calculamos la energia potencial
[Vx, S] = energia_potencial(bs, sigma, me, omega_x, x, w);
Vx = Vx(2:N_base-1,2:N_base-1);
Sx = S(2:N_base-1,2:N_base-1);

[Vy, S] = energia_potencial(bs, sigma, me, omega_y, x, w);
Vy = Vy(2:N_base-1,2:N_base-1);
Sy = S(2:N_base-1,2:N_base-1);

[Vz, S] = energia_potencial(bs, sigma, me, omega_y, x, w);
Vz = Vz(2:N_base-1,2:N_base-1);
Sz = S(2:N_base-1,2:N_base-1);

T = kron(Tx, kron(Tx, Tx));
V = kron(Vx, kron(Vy, Vz));
S = kron(Sx, kron(Sy, Sz));

estados = [];

for V0 = V0_vec
  V0 = V0/eV;

  %% Armo el hamitoniano
  H = zeros(N_base-2);
  H = T(:,:) - V0*V(:,:);

  [auvec, e] = eigs(H, S, N_base-2);
  [e, perm] = sort(diag(e));
  e = eV*e;

  auval = [V0, e'];

  auvec = auvec(:, perm);

  estados = [estados; V0, e(1), e(2), e(3), e(4), e(5), e(6)];

  save('-ascii', '-append', name, 'auval')

  % %% Ahora voy a graficar los autoestados
  % c = zeros(N_base);
  % c(2:N_base-1,2:N_base-1) = transpose(auvec);
  % bs = bspeval(kord-1, c, knots, x);
  %
  % norma = [0; diag(transpose(auvec)*auvec); 0];
  %
  % bs = bs./repmat(sqrt(norma), [1 size(bs, 2)]);
  % bs = bs.*bs;%/repmat(x, [size(bs, 1)]);
  %
  % autovectores = [x; bs]';
  % save('-ascii', 'autovectores.dat', 'autovectores')
end
%
% f = figure('visible', 'on'); hold on;
% set(gca, 'linewidth', 2, 'fontsize', 20); % grosor de la lineas de los ejes y el tamaño de la letra
% set(gca,'ticklength', 2.5*get(gca,'ticklength')) %largo de los ticksmarks
% xlabel('V0 [eV]','fontsize',20);
% ylabel('E [eV]')
%
% plot(estados(:,1), estados(:,2), '-', 'LineWidth', 3,...
%      estados(:,1), estados(:,3), '-', 'LineWidth', 3,...
%      estados(:,1), estados(:,4), '-', 'LineWidth', 3,...
%      estados(:,1), estados(:,5), '-', 'LineWidth', 3,...
%      estados(:,1), estados(:,6), '-', 'LineWidth', 3)
%
% print("1e-E_vs_V0.png")
%
% print -depsc ./E_vs_V0.eps

elapsed_time = toc()
