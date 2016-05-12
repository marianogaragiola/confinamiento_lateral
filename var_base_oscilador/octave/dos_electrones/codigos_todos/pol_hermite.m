function H = pol_hermite(N, x, omega)

  num_filas = length(x);

  H = zeros(num_filas, N+1);

  H(:,1) = sqrt(sqrt(omega/pi)); % Pol hermite de orden cero

  H(:,2) = sqrt(2*sqrt(omega/pi))*sqrt(omega)*x(:); % Pol hermite de orden uno

  % Calculo de los n-2 polinomios restantes
  for k = 3:N+1
    H(:,k) = sqrt(2/(k-1))*sqrt(omega)*(x(:).*H(:,k-1)) - sqrt((k-2)/(k-1))*H(:,k-2);
  end

end
