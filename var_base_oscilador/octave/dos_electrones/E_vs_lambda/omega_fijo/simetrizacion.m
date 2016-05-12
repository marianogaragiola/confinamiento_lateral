function H_sim = simetrizacion(H);

  N = sqrt(size(H, 1));

  [idx, idnm] = indices(N);

  N_dim = N*(N+1)/2;

  H_sim = zeros(N_dim);

  for ind1 = 1:N_dim

    n1 = idx(ind1, 2); m1 = idx(ind1, 3);

    for ind2 = 1:ind1

      n2 = idx(ind2, 2); m2 = idx(ind2, 3);

      if(m1==n1 && m2==n2)
        i1 = idnm(m1+N*(n1-1));
        j1 = idnm(m2+N*(n2-1));
        H_sim(ind1,ind2) = H(i1,j1);
      elseif(m1==n1 && m2!=n2)
        i1 = idnm(m1+N*(n1-1));
        j1 = idnm(m2+N*(n2-1)); j2 = idnm(n2+N*(m2-1));
        H_sim(ind1,ind2) = 0.5*sqrt(2)*(H(i1,j1) + H(i1,j2));
      elseif(m1!=n1 && m2==n2)
        i1 = idnm(m1+N*(n1-1)); i2 = idnm(n1+N*(m1-1));
        j1 = idnm(m2+N*(n2-1));
        H_sim(ind1,ind2) = 0.5*sqrt(2)*(H(i1,j1) + H(i1,j1));
      elseif(m1!=n1 && m2!=n2)
        i1 = idnm(m1+N*(n1-1)); i2 = idnm(n1+N*(m1-1));
        j1 = idnm(m2+N*(n2-1)); j2 = idnm(n2+N*(m2-1));
        H_sim(ind1,ind2) = 0.5*(H(i1,j1) + H(i1,j2) + H(i2,j1) + H(i2,j2));
      end

      H_sim(ind2,ind1) = H_sim(ind1,ind2);

    end

  end

  % ind1 = 1;
  % for n1 = 1, N
  %   for m1 = 1, n1
  %
  %     ind2 = 1;
  %     for n2 = 1, N
  %       for m2 = 1, n2
  %
  %         if(m1==n1 .and. m2==n2) then
  %
  %           H(ind1,ind2) = 2*H_e(n1,n2)*kronecker(n1,n2) + lambda*V_ef(n1,n1,n2,n2);
  %
  %         elseif(m1==n1 .and. m2.ne.n2) then
  %
  %           H(ind1,ind2) = sqrt2*(2*(H_e(n1,n2)*kronecker(n1,m2)+H_e(n1,m2)*kronecker(n1,n2)) +&
  %                        & lambda*(V_ef(n1,n1,n2,m2) + V_ef(n1,n1,m2,n2)));
  %
  %         elseif(m1.ne.n1 .and. m2==n2) then
  %
  %           H(ind1,ind2) = sqrt2*(2*(H_e(n1,n2)*kronecker(m1,n2)+H_e(m1,n2)*kronecker(n1,n2)) +&
  %                        & lambda*(V_ef(n1,m1,n2,n2) + V_ef(m1,n1,n2,n2)));
  %
  %         elseif(m1.ne.n1 .and. m2.ne.n2) then
  %
  %           H(ind1,ind2) = H_e(m1,m2)*kronecker(n1,n2) + H_e(n1,n2)*kronecker(m1,m2) + &
  %                        & H_e(m1,n2)*kronecker(n1,m2) + H_e(n1,m2)*kronecker(m1,n2) + &
  %                        & 0.5*lambda*(V_ef(n1,m1,n2,m2) + V_ef(n1,m1,m2,n2) + &
  %                                      &  V_ef(m1,n1,n2,m2) + V_ef(m1,n1,n2,m2) );
  %
  %         end
  %
  %         ind2 = ind2 + 1;
  %
  %       end
  %     end
  %
  %     ind1 = ind1 + 1;
  %   end
  % end

end
