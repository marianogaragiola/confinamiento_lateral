## Defino funcion que calcula los knots
import numpy as _np

def knots_sequence(grado, type, N_intervalos, beta, a, b):
	N_puntos = N_intervalos + 1;
	if type=='uniform':
		knots = _np.hstack((_np.hstack((_np.tile(a, grado), _np.linspace(a, b, N_puntos))), _np.tile(b, grado)));
	elif type=='exp':
		x = _np.linspace(a, b, N_puntos);
		y = _np.sign(x)*_np.exp(beta*_np.abs(x));
		y = (x[N_puntos-1]-x[0])/(y[N_puntos-1]-y[0])*y;
		knots = _np.hstack((_np.hstack((_np.tile(a, grado), y)), _np.tile(b, grado)));

	return knots
