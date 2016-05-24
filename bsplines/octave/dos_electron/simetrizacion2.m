function M = simetrizacion2(N)

  [idx, idnm, idsim] = indices(N);

  N_dim = N*N;

  for ind1 = 1:N_dim
    n1 = idx(ind1, 2); m1 = idx(ind1, 3);

    for ind2 = 1:N_dim
      n2 = idx(ind2, 2); m2 = idx(ind2, 3);

      M(ind1,ind2) = sqrt(0.5)*(KronD(n1,n2)*KronD(m1,m2) + KronD(n1,m2)*KronD(m1,n2))*(1-KronD(n2,m2)) + KronD(n1,n2)*KronD(m1,n2)*KronD(n2,m2);

    end
  end

end
