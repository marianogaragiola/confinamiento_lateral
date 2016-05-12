function H_sim = simetrizacion(H);

  N = sqrt(size(H, 1));

  [idx, idnm, idsim] = indices(N);

  N_dim = N*(N+1)/2;

  H_sim = zeros(N_dim);

  for ind1 = 1:N_dim

    n1 = idsim(ind1, 2); m1 = idsim(ind1, 3);

    for ind2 = 1:ind1

      n2 = idsim(ind2, 2); m2 = idsim(ind2, 3);

      % if(m1==n1 && m2==n2)
      %   i1 = idnm(m1+N*(n1-1));
      %   j1 = idnm(m2+N*(n2-1));
      %   H_sim(ind1,ind2) = H(i1,j1);
      % elseif(m1==n1 && m2~=n2)
      %   i1 = idnm(m1+N*(n1-1));
      %   j1 = idnm(m2+N*(n2-1)); j2 = idnm(n2+N*(m2-1));
      %   H_sim(ind1,ind2) = 0.5*sqrt(2)*(H(i1,j1) + H(i1,j2));
      % elseif(m1~=n1 && m2==n2)
      %   i1 = idnm(m1+N*(n1-1)); i2 = idnm(n1+N*(m1-1));
      %   j1 = idnm(m2+N*(n2-1));
      %   H_sim(ind1,ind2) = 0.5*sqrt(2)*(H(i1,j1) + H(i1,j1));
      % else %if(m1~=n1 && m2~=n2)
      %   i1 = idnm(m1+N*(n1-1)); i2 = idnm(n1+N*(m1-1));
      %   j1 = idnm(m2+N*(n2-1)); j2 = idnm(n2+N*(m2-1));
      %   H_sim(ind1,ind2) = 0.5*(H(i1,j1) + H(i1,j2) + H(i2,j1) + H(i2,j2));
      % end

      if(m1==n1 && m2==n2)
        i1 = idnm(n1,n1);
        j1 = idnm(n2,n2);
        H_sim(ind1,ind2) = H(i1,j1);
      elseif(m1==n1 && m2~=n2)
        i1 = idnm(n1,n1);
        j1 = idnm(n2,m2); j2 = idnm(m2,n2);
        H_sim(ind1,ind2) = sqrt(0.5)*(H(i1,j1) + H(i1,j2));
      elseif(m1~=n1 && m2==n2)
        i1 = idnm(n1,m1); i2 = idnm(m1,n1);
        j1 = idnm(n2,n2);
        H_sim(ind1,ind2) = sqrt(0.5)*(H(i1,j1) + H(i2,j1));
      elseif(m1~=n1 && m2~=n2)
        i1 = idnm(n1,m1); i2 = idnm(m1,n1);
        j1 = idnm(n2,m2); j2 = idnm(m2,n2);
        H_sim(ind1,ind2) = 0.5*(H(i1,j1) + H(i1,j2) + H(i2,j1) + H(i2,j2));
      end

      H_sim(ind2,ind1) = H_sim(ind1,ind2);

    end

  end
end
