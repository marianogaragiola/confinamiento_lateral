clear all;
close all;

% algunas constantes utiles
alpha = 658.4092645439; a0 = 0.0529177210; eV = 27.21138564;

% Parametros para el problema
N_cuad = 100;
N_base = 49;

sigma   = 20; % en nm
me      = 0.063; % masa efectiva del electron
B_campo = 80; % campo magnetico en tesla
omega   = 0.02; % parametro variacional no lineal
V0      = 0.05; % Profundidad del potencial en eV

lambda_vec = linspace(0, 1, 10);

% Parametros dependientes
b = 0.5*(1/sigma)^2;
l = sqrt(2*alpha/B_campo);

% empiezo con los calculos de las matrices
[x, w] = GaussHermite_2(N_cuad);

Id = eye(N_base+1);



file = fopen('./resultados/auval.dat', 'w');
fprintf(file, '# N_base = %i , N_cuad = %i \n', N_base, N_cuad);
fprintf(file, '# sigma = %f \n', sigma);
fprintf(file, '# me = %f \n', me);
fprintf(file, '# B_campo = %f \n', B_campo);
fprintf(file, '# omega = %f \n', omega);
fprintf(file, '# V0 = %f \n', V0);
fprintf(file, '# autovalores problema \n');
fclose(file);


H = Hamiltoniano(N_base, me, omega, sigma, V0, B_campo, x, w);

H = kron(H, Id) + kron(Id, H);

V_ef = Interaccion(N_base, omega, l, x, w);

for lambda = lambda_vec

  H = H + a0*lambda*V_ef;

  H = 0.5*(H+transpose(H));

  H_sim = simetrizacion(H);

  e = eig(H_sim);
  e = sort(e);
  e(:) = eV*e(:);

  resultado = [lambda, e(1:50)'];

  save('-ascii', '-append', './resultados/auval.dat', 'resultado')

end
