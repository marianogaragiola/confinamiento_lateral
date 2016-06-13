
clear all
close all
tic();

% algunas constantes utiles
a0 = 0.0529177210; eV = 27.21138564;

%% Parametros fisicos del problema
sigma   = 10; % en nm
me      = 0.063; % masa efectiva del electron

V0_vec = linspace(0, 0.03, 10);

% Parametros dependientes
b = 0.5*(1/sigma)^2;

%% Parametros para los B-splines
num_intervalos   = 12; % numero de subintervalos
num_break_points = num_intervalos + 1; % numero de break points contando los extremos
kord             = 5; % orden de los bsplines, el grado de los polinomios es kord-1
regularity       = kord-2; % repeticion de los knots, simplemente degenerados por eso kord-2
xMax             =  150; % extremo superior del intervalo
xMin             = -150; % extremo inferior del intervalo
N_cuad           = 100;
N_base = num_intervalos + kord - 1;

name = sprintf('./resultados/1e-E_vs_V0-sigma%2.0fnm.dat', sigma);

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
sigma = sigma/a0; xMax = xMax/a0; xMin = xMin/a0;

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
w = 0.5*(xMax - xMin)*w;

%% calculo los bsplines en el vector x
bs = bspeval(kord-1, c, knots, x);

%% para calcular las derivadas
[dc, dknots] = bspderiv(kord-1, c, knots);

dbs = bspeval(kord-2, dc, dknots, x);

%% Ahora paso a calcular las matrices del problema
Tx = energia_cinetica(dbs, me, w);
Tx = Tx(2:N_base-1,2:N_base-1);

%% Calculamos la energia potencial
[Vx, Sx] = energia_potencial(bs, sigma, x, w);
Vx = Vx(2:N_base-1,2:N_base-1);
Sx = Sx(2:N_base-1,2:N_base-1);

T = kron(Tx, kron(Sx, Sx)) + kron(Sx, kron(Tx, Sx)) + kron(Sx, kron(Sx, Tx));
% V = kron(Vx, kron(Sx, Sx)) + kron(Sx, kron(Vx, Sx)) + kron(Sx, kron(Sx, Vx));
V = kron(Vx, kron(Vx, Vx));
S = kron(Sx, kron(Sx, Sx));

N_dim = (N_base-2)**3;

estados = [];

for V0 = V0_vec
  V0 = V0/eV;

  %% Armo el hamitoniano
  H = zeros(N_dim);
  H = T(:,:) - V0*V(:,:);
  H = 0.5*(H + transpose(H));
  return

  [auvec, e] = eig(H, S);

  [e, perm] = sort(diag(e));
  e = eV*e;

  % auval = [eV*V0, e'];

  % auvec = auvec(:, perm);

  % estados = [estados; eV*V0, e(1), e(2), e(3), e(4), e(5), e(6)];

  % save('-ascii', '-append', name, 'auval')

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

% f = figure('visible', 'on'); hold on;
% set(gca, 'linewidth', 2, 'fontsize', 20); % grosor de la lineas de los ejes y el tamaño de la letra
% set(gca,'ticklength', 2.5*get(gca,'ticklength')) %largo de los ticksmarks
% xlabel('V0 [eV]','fontsize',20);
% ylabel('E [eV]')
%
% plot(estados(:,1), estados(:,2), '-', 'LineWidth', 3,...
%     estados(:,1), estados(:,3), '-', 'LineWidth', 3,...
%     estados(:,1), estados(:,4), '-', 'LineWidth', 3,...
%     estados(:,1), estados(:,5), '-', 'LineWidth', 3,...
%     estados(:,1), estados(:,6), '-', 'LineWidth', 3)
%
% print -depsc -color 1e-E_vs_V0.eps

elapsed_time = toc()
