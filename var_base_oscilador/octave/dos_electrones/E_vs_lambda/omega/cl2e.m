clear all;
close all;

% algunas constantes utiles
alpha = 658.4092645439; a0 = 0.0529177210; eV = 27.21138564;

% Parametros para el problema
N_cuad = 100;
N_base = 39;

sigma   = 20; % en nm
me      = 0.063; % masa efectiva del electron
B_campo = 80; % campo magnetico en tesla
omega   = 0.02; % parametro variacional no lineal
V0      = 0.05; % Profundidad del potencial en eV
N_omega = 10;

lambda_vec = linspace(0, 1, 10);
omega_vec = linspace(1e-4, 0.1, N_omega);

% Parametros dependientes
b = 0.5*(1/sigma)^2;
l = sqrt(2*alpha/B_campo);

% empiezo con los calculos de las matrices
[x, w] = GaussHermite_2(N_cuad);

Id = eye(N_base+1);

for ind = 1:N_omega

  omega = omega_vec(ind);

  name = sprintf('./resultados/auval_vs_lambda-B%2.0f-omega%d.dat', B_campo, ind);

  file = fopen(name, 'w');
  fprintf(file, '# N_base = %i , N_cuad = %i \n', N_base, N_cuad);
  fprintf(file, '# sigma = %f nm\n', sigma);
  fprintf(file, '# me = %f UA\n', me);
  fprintf(file, '# B_campo = %f T\n', B_campo);
  fprintf(file, '# omega = %f \n', omega);
  fprintf(file, '# V0 = %f eV\n', V0);
  fprintf(file, '# autovalores problema \n');
  fclose(file);

  H = Hamiltoniano(N_base, me, omega, sigma, V0, B_campo, x, w);

  % H = kron(H, Id) + kron(Id, H);

  V_ef = Interaccion(N_base, omega, l, x, w);

  for lambda = lambda_vec

    % H = H + lambda*V_ef;

    % H = 0.5*(H+transpose(H));

    % H_sim = simetrizacion(H);
    H_sim = simetrizacion2(H, V_ef, lambda);

    e = eig(H_sim);
    e = sort(e);
    e(:) = eV*e(:);

    resultado = [lambda, e(1:50)'];

    save('-ascii', '-append', name, 'resultado')

  end

end
