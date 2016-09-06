# Codigo que calcula los autovalores y autovectores del hamiltoniano
# efectio 1D para un electron.
#
# Terminar la descripcion
#
#######################################################################
#######################################################################
import numpy as np
import matplotlib.pyplot as plt
from scipy import linalg as LA
from bspline import Bspline
from sys import exit

## Defino funcion que calcula los knots
def knots_sequence(grado, type, N_intervalos, beta, a, b):
	N_puntos = N_intervalos + 1;
	if type=='uniform':
		knots = np.hstack((np.hstack((np.tile(a, grado), np.linspace(a, b, N_puntos))), np.tile(b, grado)));
	elif type=='exp':
		x = np.linspace(a, b, N_puntos);
		y = np.sign(x)*np.exp(beta*np.abs(x));
		y = (x[N_puntos-1]-x[0])/(y[N_puntos-1]-y[0])*y;
		knots = np.hstack((np.hstack((np.tile(a, grado), y)), np.tile(b, grado)));

	return knots
#############################

## Defino el potencial del pozo
def V_potencial_pozo(sigma, x):
	U = np.exp(-0.5*x**2/sigma**2);
	return U

#######################################
#######################################
## Aca empieza el main del programa ###
#######################################
#######################################

## Parametros utiles
a0 = 0.0529177210; eV = 27.21138564; c_light = 137.035999074492; ua_to_T = 1.72e3;
alpha = 658.4092645439;

## Parametros fisicos del problema
me = 0.063; ## Masa de la particula
sigma = 2.0; ## Ancho del pozo gaussiano en nm
v1 = 0.05; ## Parametro del pozo en eV/nm^2
v2 = 0.45; ## Intensidad de campo magnetico en teslas

bcampo_vec = np.linspace(1.0, 1000.0, 99);

## Intervalo de integracion en la coordenada z
Zmax = 100.0;
Zmin = -Zmax;

N_intervalos_z = 50;
N_cuad = 100;
grado = 4;
kord = grado + 1;
beta = 0.0065; ## Cte de decaimiento para los knots en la distribucion exponencial
N_splines_z = N_intervalos_z + grado;
N_base_z = N_splines_z - 2;

N_dim = N_base_z; # Tamano de las matrices

nev = N_dim;

if nev > N_dim:
	print "nev > N_dim, cambiar el valor de nev"
	exit(0)

archivo = "./resultados/1e-E_vs_B-v1_%6.4feVsnm2-v2_%6.4feV-sigma_%6.4f.dat" % (v1, v2, sigma)

f = open(archivo, 'w')
f.write("# Intervalo de integracion en z [{0}, {1}]\n".format(Zmin, Zmax))
f.write("# Grado de los B-splines {0}\n".format(grado))
f.write("# Num de intervalos {0} y tamano de base {1} en z\n".format(N_intervalos_z, N_base_z))
f.write("# Dimension total del espacio N_dim = N_base_z = {0}\n".format(N_dim))
f.write("# Cte de separacion de knots en dist exp beta = {0}\n".format(beta))
f.write("# Orden de la cuadratura N_cuad = {0}\n".format(N_cuad))
f.write("# Masa de la particula me = {0} en UA\n".format(me))
f.write("# Ancho del pozo gaussiano sigma = {0} en nm\n".format(sigma))
f.write("# Parametro del pozo v1 = {0} en eV/nm^2\n".format(v1))
f.write("# Parametro del pozo v2 = {0} en eV\n".format(v2))
f.write("# Autovalores calculados\n")

## Paso a unidades atomicas
Zmax = Zmax/a0; Zmin = Zmin/a0;
sigma = sigma/a0; beta = beta*a0;
v1 = v1*a0**2/eV; v2 = v2/eV;

## Vector de knots para definir los B-splines, distribucion uniforme
knots_z = knots_sequence(grado, 'uniform', N_intervalos_z, beta, Zmin, Zmax);

## Pesos y nodos para la cuadratura de Gauss-Hermite
x, w = np.polynomial.legendre.leggauss(N_cuad);

## Escribo bien los nodos y los pesos para las integrales en z
z_nodos = np.array([]);
wz_pesos = np.array([]);
for i in range(N_intervalos_z+1):
	aux_x = 0.5*(knots_z[i+grado+1]-knots_z[i+grado])*x + 0.5*(knots_z[i+grado+1]+knots_z[i+grado]);
	aux_w = 0.5*(knots_z[i+grado+1]-knots_z[i+grado])*w;

	z_nodos = np.hstack((z_nodos, aux_x));
	wz_pesos = np.hstack((wz_pesos, aux_w));

wz_pesos = np.tile(wz_pesos, (N_splines_z, 1));
# z_nodos2 = np.tile(z_nodos, (N_splines_z, 1));

## B-splines en la coordenada z
basis = Bspline(knots_z, grado);
## Calculo la matriz con los bsplines y sus derivadas en z
bsz  = [basis._Bspline__basis(i, basis.p) for i in z_nodos]; # evaluo los bsplines
dbsz = [basis.d(i) for i in z_nodos];                        # evaluo las derivadas de los bsplines

## Matriz de solapamiento en z
Sz = np.dot(np.transpose(bsz), (np.transpose(wz_pesos)*bsz));
Sz = np.array([[Sz[i][j] for i in range(1, N_splines_z-1)] for j in range(1, N_splines_z-1)]);

## Matriz de energia cinetica en z
Tz = 0.5/me*np.dot(np.transpose(dbsz), (np.transpose(wz_pesos)*dbsz));
Tz = np.array([[Tz[i][j] for i in range(1, N_splines_z-1)] for j in range(1, N_splines_z-1)]);

## Matriz del potencial de confinamiento
expz = V_potencial_pozo(sigma, z_nodos);

Uz = z_nodos**2*expz;

expz = np.tile(expz, (N_splines_z, 1));
Uz = np.tile(Uz, (N_splines_z, 1));

V1z = np.dot(np.transpose(bsz), (np.transpose(expz)*np.transpose(wz_pesos)*bsz));
V1z = np.array([[V1z[i][j] for i in range(1, N_splines_z-1)] for j in range(1, N_splines_z-1)]);

V2z = np.dot(np.transpose(bsz), (np.transpose(Uz)*np.transpose(wz_pesos)*bsz));
V2z = np.array([[V2z[i][j] for i in range(1, N_splines_z-1)] for j in range(1, N_splines_z-1)]);

## Calculo los autovalores y autovectores
auval = np.zeros(N_dim);
for bcampo in bcampo_vec:

	l_campo = np.sqrt(2.0*alpha/bcampo)/a0;
	b = 0.5/sigma**2;

	cte1 = l_campo**2/(1. + l_campo**2*b);
	cte2 = 1./(1. + l_campo**2*b);

	v1_2 = v1/(1.0 + l_campo**2*b)
	v2_2 = (v1*l_campo**2 - v2)/(1.0 + l_campo**2*b)

	# El hamiltoniano es
	# Hz = Tz + (v1*cte1-v2*cte2)*V1z + v1*cte2*V2z
	Hz = Tz + v1_2*V2z + v2_2*V1z

	e, auvec = LA.eigh(Hz, Sz);
	e = eV*e;

	auval = np.vstack((auval, e));

	f.write("{:13.9e}   ".format(bcampo))
	for i in range(nev):
		f.write("{:13.9e}   ".format(e[i]))

	f.write("\n")

f.close()
