function V_ef = Interaccion(N_base, omega, l, x, w)

  a0 = 0.0529177210;

  Hermite = pol_hermite(N_base, x, 1);

  % E = kron(Hermite, Hermite);
  %
  % Her_2 = E([1:length(x)+1:(length(x))^2],:);
  % Her_2 = repmat(w, [1 size(Her_2, 2)]).*Her_2;
  %
  % x_2 = repmat(x, [1 length(x)]);
  % y = (x_2-transpose(x_2))/(sqrt(2*omega)*l);
  %
  % pot = sqrt(0.5*pi)/l*erfcx(abs(y));
  % % pot = (repmat(w, [1 length(w)]).*repmat(transpose(w), [length(w) 1])).*pot;
  %
  % V_ef = transpose(Her_2)*(pot*Her_2);
  % % V_ef = 0.5*(V_ef + transpose(V_ef));


  [idx, idnm, idsim] = indices(N_base+1);

  x_2 = repmat(x, [1 length(x)]);
  y = (x_2-transpose(x_2))/(sqrt(2*omega)*l);

  pot = sqrt(0.5*pi)/l*erfcx(abs(y))*a0;

  N_dim = (N_base+1)^2;
  V_ef = zeros(N_dim);
  for ind1 = 1:N_dim

    n1 = idx(ind1, 2); m1 = idx(ind1, 3);

    for ind2 = 1:ind1

      n2 = idx(ind2, 2); m2 = idx(ind2, 3);

      if(mod(n1+n2,2)==0 && mod(m1+m2,2)==0)
        phi_n1n2 = w(:).*Hermite(:,n1).*Hermite(:,n2);
        phi_m1m2 = w(:).*Hermite(:,m1).*Hermite(:,m2);

        V_ef(ind1,ind2) = transpose(phi_n1n2)*(pot*phi_m1m2);
      end

      V_ef(ind2,ind1) = V_ef(ind1,ind2);

    end
  end

end
