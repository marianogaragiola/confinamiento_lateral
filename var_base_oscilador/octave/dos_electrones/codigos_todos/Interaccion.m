function V_ef = Interaccion(N_base, omega, l, x, w)

  Hermite = pol_hermite(N_base, x, 1);

  E = kron(Hermite, Hermite);

  Her_2 = E([1:length(x)+1:(length(x))^2],:);
  Her_2 = repmat(w, [1 size(Her_2, 2)]).*Her_2;

  x_2 = repmat(x, [1 length(x)]);
  y = (x_2-transpose(x_2))/(sqrt(2*omega)*l);

  pot = sqrt(0.5*pi)/l*erfcx(abs(y));
  % pot = (repmat(w, [1 length(w)]).*repmat(transpose(w), [length(w) 1])).*pot;

  V_ef = transpose(Her_2)*(pot*Her_2);
  V_ef = 0.5*(V_ef + transpose(V_ef));

end
