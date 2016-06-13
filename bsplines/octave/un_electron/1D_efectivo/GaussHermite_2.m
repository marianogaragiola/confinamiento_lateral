function [x, w] = GaussHermite_2(n)

i = 1:n-1;
a = sqrt(i/2);
CM = diag(a,1) + diag(a,-1);

[V L] = eig(CM);
[x ind] = sort(diag(L));
V = V(:,ind)';
w = sqrt(pi) * V(:,1).^2;

end
