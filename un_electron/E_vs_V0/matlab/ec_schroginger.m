clear all;
close all;
format long;

eV = 27.21138564;
me = 0.036;
sigma = 15;
N = 640;
z_i = -80.0; z_f = 80.0;
z = linspace(z_i, z_f, N);
z = z(:);
V0_i = 0.0; V0_f = 10; % en eV
V = linspace(V0_i, V0_f, 150);

Dz = (z_f-z_i)/(N-1);

% Armo la energia potencial
f = -exp(-z.*z/(2.0*sigma^2));
potencial = diag(f);

% Armo la matrix de la energia cinetica
T = diag(repmat(-2, [N-2,1])) + diag(repmat(1, [N-3,1]), 1) + diag(repmat(1, [N-3,1]), -1);
cinetica = -T/Dz^2;

for V0 = V
  % Armo el hamiltiano
  H = eV*cinetica + V0*potencial(2:N-1,2:N-1);

  % Acordarse que el vector es u = (u1, u2, ....., uN), pero las condiciones de contorno
  % son u1 = un = 0, por lo tanto solo calculo u = (u2, u3, ..., uN-1)

  % Ahora calculo los autovalores y autovectores de M, que deberian ser las soluciones a la ecuacion
  % de onda de una cuerda de largo L

  [auvec_aux, E] = eig(H);

  energias = diag(E)';
  energias = [V0, energias];

  auvec = zeros(N);

  auvec(2:N-1,2:N-1) = auvec_aux;

  autoestados = auvec(:,1:10);

  autoestados = [z, autoestados];

  save('-ascii', '-append', './auval.dat', 'energias');

end
