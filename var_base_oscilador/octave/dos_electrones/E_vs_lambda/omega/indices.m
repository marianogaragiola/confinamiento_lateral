% function [idx, idnm] = indices(N)
%
%   M = N*(N+1)/2;
%
%   idx = [(1:1:M)'];
%
%   aux = [];
%   for n = 1:N
%     aux = [aux; repmat(n, [n 1]) (1:1:n)'];
%   end
%
%   idx = [idx aux];
%
%   idnm = [(1:1:N^2)];
%
% end

function [idx, idnm, idsim] = indices(N)

  M = N*(N+1)/2;

  idx = [(1:1:N*N)'];
  aux = [];
  idnm = [];
  for n = 1:N
    aux = [aux; repmat(n, [N 1])];
    idnm = [idnm; ((n-1)*N+1:1:n*N)];
  end

  idx = [idx aux repmat(transpose((1:1:N)), [N 1])];

  idsim = [(1:1:M)'];

  aux = [];
  for n = 1:N
    aux = [aux; repmat(n, [n 1]) (1:1:n)'];
  end

  idsim = [idsim aux];

end
