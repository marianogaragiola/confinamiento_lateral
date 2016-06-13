function T = energia_cinetica(dbs, me, w);

  T = 0.5/me*(dbs.*repmat(w, [size(dbs,1) 1]))*transpose(dbs);

end
