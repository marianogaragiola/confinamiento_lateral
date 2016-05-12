function [V, Id] = energia_potencial(N_base, omega, sigma, x, w)

  Her = pol_hermite(N_base, x, 1);

  V = zeros(N_base+1);

  potencial = zeros([length(x) 1]);
  potencial(:) = exp(-0.5/omega*(x(:)/sigma).^2);

  potencial(:) = w(:).*potencial(:);

  V = transpose(Her)*(repmat(potencial, [1, N_base+1]).*Her);

  Id = transpose(Her)*(repmat(w, [1, N_base+1]).*Her);
end
