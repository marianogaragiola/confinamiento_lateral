function [idx, idnm] = indices(N)

  M = N*(N+1)/2;

  idx = [(1:1:M)'];

  aux = [];
  for n = 1:N
    aux = [aux; repmat(n, [n 1]) (1:1:n)'];
  end

  idx = [idx aux];

  idnm = [(1:1:N^2)];

end
