clear all
close all

n = 3; m = 4;
H = rand(n,m);

E = kron(H, H);
F = E([1:n+1:n*n],:);
