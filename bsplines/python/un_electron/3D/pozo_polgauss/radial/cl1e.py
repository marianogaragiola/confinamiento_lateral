# Codigo en python que calcula los autovalores de una
# particula en un pozo guassiano de simetrica esferica
# al cual se le aplica un campo magnetico.
# Se calcula los elementos de matriz del hamiltoniano
# en coordenadas cilindricas.
# El hamiltoniano es
#
# H = -hbar^2/(2*me) (1/r*d/dr(r*d/dr)+d^2/dz^2) + V(r,z) +
#     + 0.5*me*omega^2*r^2 + hbar^2m^2/(2*me*r^2) + hbar*omega*m
#
# donde m es la componente en el eje z del momento angular.
#
# V(r,z) = -v0*exp(-(r^2+z^2)/(2*sigma))
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
def V_potencial_pozo(v1, v2, sigma, x):
	U = (v1*x**2-v2)*np.exp(-0.5*x**2/sigma**2);
	return U

#############################

## Energia cinetica del momento angular
def L_angular(me, l_angular, x):
	L = 0.5*l_angular*(l_angular+1)/(me*x**2);
	return L

#######################################
#######################################
## Aca empieza el main del programa ###
#######################################
#######################################

## Parametros utiles
a0 = 0.0529177210; eV = 27.21138564; c_light = 137.035999074492; ua_to_T = 1.72e3;

## Parametros fisicos del problema
me = 0.063; ## Masa de la particula
l_angular = 0.0; ## Momento angular
sigma = 2.0; ## Ancho del pozo gaussiano en nm
v1 = 0.05; ## Parametro del pozo en eV/nm^2

v2_vec = np.linspace(0.0, 1.0, 100);

## Separo las coordenadas y tomo distinta base en r y en z
Rmin = 0.0;
Rmax = 100.0;

N_intervalos_r = 30;
N_cuad = 200;
grado = 4;
kord = grado + 1;
beta = 0.0065; ## Cte de decaimiento para los knots en la distribucion exponencial
N_splines_r = N_intervalos_r + grado;
N_base_r = N_splines_r - 2;

N_dim = N_base_r; # Tamano de las matrices

archivo = "./resultados/1e-E_vs_v2-v1_%6.4fT-l_%d.dat" % (v1, l_angular)

f = open(archivo, 'w')
f.write("# Intervalo de integracion en r [{0}, {1}]\n".format(Rmin, Rmax))
f.write("# Grado de los B-splines {0}\n".format(grado))
f.write("# Num de intervalos {0} y tamano de base {1} en r\n".format(N_intervalos_r, N_base_r))
f.write("# Dimension total del espacio N_dim = N_base_r = {0}\n".format(N_dim))
f.write("# Cte de separacion de knots en dist exp beta = {0}\n".format(beta))
f.write("# Orden de la cuadratura N_cuad = {0}\n".format(N_cuad))
f.write("# Masa de la particula me = {0} en UA\n".format(me))
f.write("# Momento angular l_angular = {0}\n".format(l_angular))
f.write("# Ancho del pozo gaussiano sigma = {0} en nm\n".format(sigma))
f.write("# Parametro v1 del pozo de potencial v1 = {} en eV/nm^2\n".format(v1))
f.write("# Autovalores calculados\n")

## Paso a unidades atomicas
Rmax = Rmax/a0; Rmin = Rmin/a0;
sigma = sigma/a0; beta = beta*a0;
v1 = v1*a0**2/eV;

## Vector de knots para definir los B-splines, distribucion uniforme
knots_r = knots_sequence(grado, 'uniform', N_intervalos_r, beta, Rmin, Rmax);

## Pesos y nodos para la cuadratura de Gauss-Hermite
x, w = np.polynomial.legendre.leggauss(N_cuad);

## Escribo bien los nodos y los pesos para las integrales en r
r_nodos = np.array([]);
wr_pesos = np.array([]);
for i in range(N_intervalos_r+1):
	aux_x = 0.5*(knots_r[i+grado+1]-knots_r[i+grado])*x + 0.5*(knots_r[i+grado+1]+knots_r[i+grado]);
	aux_w = 0.5*(knots_r[i+grado+1]-knots_r[i+grado])*w;

	r_nodos = np.hstack((r_nodos, aux_x));
	wr_pesos = np.hstack((wr_pesos, aux_w));

wr_pesos = np.tile(wr_pesos, (N_splines_r, 1));

## B-splines en la coordenada r
basis = Bspline(knots_r, grado);
## Calculo la matriz con los bsplines y sus derivadas en r
bsr  = [basis._Bspline__basis(i, basis.p) for i in r_nodos]; # evaluo los bsplines
dbsr = [basis.d(i) for i in r_nodos];                        # evaluo las derivadas de los bsplines


## Matriz de solapamiento en r
Sr = np.dot(np.transpose(bsr), (np.transpose(wr_pesos)*bsr));
Sr = np.array([[Sr[i][j] for i in range(1,N_splines_r-1)] for j in range(1,N_splines_r-1)]);

## Matriz de energia cinetica en r
Tr = 0.5/me*np.dot(np.transpose(dbsr), (np.transpose(wr_pesos)*dbsr));
Tr = np.array([[Tr[i][j] for i in range(1,N_splines_r-1)] for j in range(1,N_splines_r-1)]);

## Matriz de momento angular al cuadrado
L = L_angular(me, l_angular, r_nodos);
L = np.tile(L, (N_splines_r, 1));
VL = np.dot(np.transpose(bsr), (np.transpose(L)*np.transpose(wr_pesos)*bsr));
VL = np.array([[VL[i][j] for i in range(1,N_splines_r-1)] for j in range(1,N_splines_r-1)]);

auval = np.zeros(N_dim);
for v2 in v2_vec:
	v2 = v2/eV;

	## Matriz de energia de potencial del pozo de potencial
	# Primero en la variable r
	U = V_potencial_pozo(v1, v2, sigma, r_nodos);
	U = np.tile(U, (N_splines_r, 1));
	Vp_r = np.dot(np.transpose(bsr), (np.transpose(U)*np.transpose(wr_pesos)*bsr));
	Vp_r = np.array([[Vp_r[i][j] for i in range(1,N_splines_r-1)] for j in range(1,N_splines_r-1)]);
	
	## El hamiltoniano en la coordenada r es
	Hr = Tr + VL + Vp_r;

	e, auvec = LA.eigh(Hr, Sr);
	e = eV*e;

	auval = np.vstack((auval, e));

	f.write("{:13.9e}   ".format(eV*v2))
	for i in range(N_dim):
		f.write("{:13.9e}   ".format(e[i]))

	f.write("\n")

auval = np.array([[auval[i][j] for i in range(1,np.size(v2_vec)+1)] for j in range(N_dim)]);

for i in range(10):
	estado = np.zeros(np.size(v2_vec));
	for j in range(np.size(v2_vec)):
 		estado[j] = auval[i][j];

 	plt.plot(v2_vec, estado, '-')

plt.show()
# 
# auval = np.vstack((np.transpose(v0_vec), auval));
#
# np.savetxt(archivo, np.transpose(auval))
f.close()
