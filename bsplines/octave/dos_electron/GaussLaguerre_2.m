
function [x, w] = GaussLaguerre_2(n, alpha)

i = 1:n;
a = (2*i-1) + alpha;
b = sqrt( i(1:n-1) .* ((1:n-1) + alpha) );
CM = diag(a) + diag(b,1) + diag(b,-1);

[V L] = eig(CM);
[x ind] = sort(diag(L));
V = V(:,ind)';
w = gamma(alpha+1) .* V(:,1).^2;

end
