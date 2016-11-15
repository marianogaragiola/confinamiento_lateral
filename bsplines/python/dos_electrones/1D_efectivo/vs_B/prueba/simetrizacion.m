function [H_sim, S_sim] = simetrizacion(H, S, V_ef, lambda);

  N = size(H, 1);

  [idx, idnm, idsim] = indices(N);

  N_dim = N*(N+1)/2;

  H_sim = zeros(N_dim);
  S_sim = zeros(N_dim);

  for ind1 = 1:N_dim

    n1 = idsim(ind1, 2); m1 = idsim(ind1, 3);

    for ind2 = 1:ind1

      n2 = idsim(ind2, 2); m2 = idsim(ind2, 3);

      if(m1==n1 && m2==n2)
        i1 = idnm(n1,n1);
        j1 = idnm(n2,n2);
        H_sim(ind1,ind2) = 2*H(n1,n2)*S(n1,n2) + lambda*V_ef(i1,j1);
        S_sim(ind1,ind2) = S(n1,n2)*S(n1,n2);
      elseif(m1==n1 && m2~=n2)
        i1 = idnm(n1,n1);
        j1 = idnm(n2,m2); j2 = idnm(m2,n2);
        H_sim(ind1,ind2) = sqrt(0.5)*(2*(H(n1,n2)*S(n1,m2) + H(n1,m2)*S(n1,n2)) ...
                         + lambda*(V_ef(i1,j1) + V_ef(i1,j2)));

        S_sim(ind1,ind2) = sqrt(0.5)*2*S(n1,n2)*S(n1,m2);
      elseif(m1~=n1 && m2==n2)
        i1 = idnm(n1,m1); i2 = idnm(m1,n1);
        j1 = idnm(n2,n2);
        H_sim(ind1,ind2) = sqrt(0.5)*(2*(H(m1,n2)*S(n1,n2) + H(n1,n2)*S(m1,n2)) ...
                         + lambda*(V_ef(i1,j1) + V_ef(i2,j1)));

        S_sim(ind1,ind2) = sqrt(0.5)*2*S(n1,n2)*S(m1,n2);
      elseif(m1~=n1 && m2~=n2)
        i1 = idnm(n1,m1); i2 = idnm(m1,n1);
        j1 = idnm(n2,m2); j2 = idnm(m2,n2);
        H_sim(ind1,ind2) = H(m1,m2)*S(n1,n2) + H(n1,n2)*S(m1,m2) ...
                         + H(m1,n2)*S(n1,m2) + H(n1,m2)*S(m1,n2) ...
                         + 0.5*lambda*(V_ef(i1,j1) + V_ef(i1,j2) + V_ef(i2,j1) + V_ef(i2,j2));

        S_sim(ind1,ind2) = S(n1,n2)*S(m1,m2) + S(n1,m2)*S(m1,n2);

      end

      H_sim(ind2,ind1) = H_sim(ind1,ind2);
      S_sim(ind2,ind1) = S_sim(ind1,ind2);

    end

  end
end
