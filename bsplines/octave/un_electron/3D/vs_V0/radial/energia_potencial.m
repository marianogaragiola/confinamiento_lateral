function [V, S] = energia_potencial(bs, sigma, me, l_ang, V0, x, w)

  %% Armo la matriz de solapamiento
  S = (bs.*repmat(w, [size(bs,1) 1]))*transpose(bs);

  pot = -V0*exp(-0.5*(x.*x)/sigma^2) + 0.5*l_ang*(l_ang + 1)./(me*x.*x);
  pot = w.*pot;

  V = (bs.*repmat(pot, [size(bs, 1) 1]))*transpose(bs);

end
