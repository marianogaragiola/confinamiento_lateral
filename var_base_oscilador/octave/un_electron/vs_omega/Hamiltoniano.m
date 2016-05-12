function H = Hamiltoniano(N_base, me, omega, sigma, V0, B_campo, x, w)

  alpha = 658.4092645439; a0 = 0.0529177210; eV = 27.21138564;

  % Parametros dependientes
  b = 0.5*(1/sigma)^2;
  l = sqrt(2*alpha/B_campo);

  % calculamos las matrices de la energia cinetica y potencial
  T = energia_cinetica(N_base, me);

  [V, Id] = energia_potencial(N_base, omega, sigma, x, w);

  H = a0^2*omega*T - V0/eV*V/(1 + l^2*b);

end
