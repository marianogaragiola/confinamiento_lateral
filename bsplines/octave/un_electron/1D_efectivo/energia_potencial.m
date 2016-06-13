function [V, S] = energia_potencial(bs, sigma, x, w)

  %% Armo la matriz de solapamiento
  S = (bs.*repmat(w, [size(bs,1) 1]))*transpose(bs);

  pot = exp(-0.5*(x.*x)/sigma^2);
  pot = w.*pot;

  V = (bs.*repmat(pot, [size(bs, 1) 1]))*transpose(bs);

end
