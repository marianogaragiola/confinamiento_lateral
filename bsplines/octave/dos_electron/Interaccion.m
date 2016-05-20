function V_ef = Interaccion(bs, x, w, l)

  N_base = size(bs, 1) - 2; % Esto porque no considero el primer y ultimo spline

  bs = transpose(bs);
  % x = transpose(x); w = transpose(w);

  [idx, idnm, idsim] = indices(N_base);

  x_2 = repmat(x, [length(x) 1]);
  y = (x_2-transpose(x_2))/(sqrt(2)*l);

  E = kron(bs, bs);

  bs_2 = E([1:length(x)+1:(length(x))^2],:);
  bs_2 = repmat(transpose(w), [1 size(bs_2, 2)]).*bs_2;

  pot = sqrt(0.5*pi)/l*erfcx(abs(y));

  V_ef = transpose(bs_2)*(pot*bs_2);
  % V_ef = 0.5*(V_ef + transpose(V_ef));

  % pot = sqrt(0.5*pi)/l*erfcx(abs(y));
  %
  % N_dim = N_base^2;
  %
  % V_ef = zeros(N_dim);
  %
  % for ind1 = 1:N_dim
  %
  %   n1 = idx(ind1, 2) + 1; m1 = idx(ind1, 3) + 1;
  %
  %   for ind2 = 1:ind1
  %
  %     n2 = idx(ind2, 2) + 1; m2 = idx(ind2, 3) + 1;
  %
  %     if(mod(n1+n2,2)==0 && mod(m1+m2,2)==0)
  %       phi_n1n2 = w(:).*bs(:,n1).*bs(:,n2);
  %       phi_m1m2 = w(:).*bs(:,m1).*bs(:,m2);
  %
  %       V_ef(ind1,ind2) = transpose(phi_n1n2)*(pot*phi_m1m2);
  %     end
  %
  %     V_ef(ind2,ind1) = V_ef(ind1,ind2);
  %
  %   end
  % end

end
