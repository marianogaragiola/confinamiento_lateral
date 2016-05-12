function T = energia_cinetica(N_dim, me)

  T = 0;
  
  diagonal = 0:1:N_dim;

  diag_up = sqrt(diagonal.*(diagonal-1));

  diag_down = sqrt((diagonal+1).*(diagonal+2));

  T = diag(diagonal+0.5) - 0.5*diag(diag_up(3:N_dim+1), 2) - 0.5*diag(diag_down(1:N_dim-1), -2);
  T = 0.5*T/me;

end
