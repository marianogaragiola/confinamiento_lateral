# Codigo en python que calcula los autovalores de una
# particula en un pozo de simetrica esferica
# al cual se le aplica un campo magnetico.
# Se calcula los elementos de matriz del hamiltoniano
# en coordenadas cilindricas.
# El hamiltoniano es
#
# H = -hbar^2/(2*me) (1/r*d/dr(r*d/dr)+d^2/dz^2) + V(r,z) +
#     + 0.5*me*omega^2*r^2 + hbar^2mz^2/(2*me*r^2) + hbar*omega*mz
#
# donde mz es la componente en el eje z del momento angular.
#
# El potencial de confinamiento para la particula es
#
#           v1 if(r<r0 and az/2<|z|<(az+bz)/2)
# V(r,z) = -v2 if(r>r0 and |z|<az/2)
#           0
#
# La matriz del hamiltoniano la calculo usando como base
# del espacio los B-splines. El problema de autovalores es un
# problema generalizado
#
#     H*c = E*S*c
#
# Guardo, los autovalores del hamiltoniano en funcion del
# campo magnetico aplicado.
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

## Defino el potencial del pozo en la coordenada rho
def V_potencial_pozo_rho(r0, x):
	U = np.piecewise(x, [x<r0, r0<=x], [1, 0]);
	return U
#############################

#######################################
#######################################
## Aca empieza el main del programa ###
#######################################
#######################################

## Parametros utiles
a0 = 0.0529177210; eV = 27.21138564; c_light = 137.035999074492; ua_to_T = 1.72e3;

## Parametros fisicos del problema
me		 = 0.041; #0.063; ## Masa de la particula
mz		 = 0.0; ## Componente z del momento angular
r0		 = 7.0; ## Ancho del pozo en rho
v0		 = 0.108844; ## Profundidad del pozo


## Separo las coordenadas y tomo distinta base en r y en z
Rmin = 0.0;
Rmax = 50.0;

N_intervalos_r = 50;
N_cuad = 500;
grado = 6;
kord = grado + 1;
beta = 0.0065; ## Cte de decaimiento para los knots en la distribucion exponencial
N_splines_r = N_intervalos_r + grado;
N_base_r = N_splines_r - 1;

N_dim = N_base_r # Tamano de las matrices

archivo1 = "./resultados/E-v0_%6.4feV-r0_%6.4fnm.dat" % (v0, r0)
archivo2 = "./resultados/densidad-v0_%6.4feV-r0_%6.4fnm.dat" % (v0, r0)

f1 = open(archivo1, 'w')
f2 = open(archivo2, 'w')
f1.write("# Intervalo de integracion en r [{0}, {1}]\n".format(Rmin, Rmax))
f1.write("# Grado de los B-splines {0}\n".format(grado))
f1.write("# Num de intervalos {0} y tamano de base {1} en r\n".format(N_intervalos_r, N_base_r))
f1.write("# Dimension total del espacio N_dim = N_base_r = {0}\n".format(N_dim))
f1.write("# Cte de separacion de knots en dist exp beta = {0}\n".format(beta))
f1.write("# Orden de la cuadratura N_cuad = {0}\n".format(N_cuad))
f1.write("# Masa de la particula me = {0} en UA\n".format(me))
f1.write("# Componente z del momento angular mz = {0}\n".format(mz))
f1.write("# Ancho del pozo en rho r_0 = {0} en nm\n".format(r0))
f1.write("# Profundidad del pozo v0 = {0} en eV\n".format(v0))
f1.write("# Autovalores calculados\n")

## Paso a unidades atomicas
Rmax = Rmax/a0; Rmin = Rmin/a0;
beta = beta*a0;
r0 = r0/a0;
v0 = v0/eV;

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

w_pesos = wr_pesos
wr_pesos = np.tile(wr_pesos, (N_splines_r-1, 1));
r_nodos2 = np.tile(r_nodos, (N_splines_r-1, 1));

## B-splines en la coordenada r
basis = Bspline(knots_r, grado);
## Calculo la matriz con los bsplines y sus derivadas en r
bsr  = [basis._Bspline__basis(i, basis.p) for i in r_nodos]; # evaluo los bsplines
dbsr = [basis.d(i) for i in r_nodos];                        # evaluo las derivadas de los bsplines

bsr = np.array(bsr)
bsr = bsr[:,0:N_splines_r-1]

dbsr = np.array(dbsr)
dbsr = dbsr[:,0:N_splines_r-1]

## Matriz de solapamiento en r
Sr = np.dot(np.transpose(bsr), (np.transpose(r_nodos2)*np.transpose(wr_pesos)*bsr));
# Sr = np.array([[Sr[i][j] for i in range(N_splines_r-1)] for j in range(N_splines_r-1)]);

## Matriz de energia cinetica en r
Tr = 0.5/me*np.dot(np.transpose(dbsr), (np.transpose(r_nodos2)*np.transpose(wr_pesos)*dbsr));
#Tr = np.array([[Tr[i][j] for i in range(N_splines_r-1)] for j in range(N_splines_r-1)]);

## Matriz de energia de potencial del pozo de potencial
Ur = V_potencial_pozo_rho(r0, r_nodos)
Ur = np.tile(Ur, (N_splines_r-1, 1))

Vr = np.dot(np.transpose(bsr), (np.transpose(r_nodos2*Ur*wr_pesos)*bsr));
# Vr = np.array([[Vr[i][j] for i in range(N_splines_r-1)] for j in range(N_splines_r-1)]);

Hr = Tr - v0*Vr

e, auvec = LA.eigh(Hr, Sr);
e = eV*e;

ge = auvec[:,0]

dge = np.dot(dbsr, ge)
ge = np.dot(bsr, ge)

# ge = r_nodos*ge**2 #/a0

file1 = "ground_state.dat"
file2 = "d_ground_state.dat"

f3 = open(file1, 'w')
f4 = open(file2, 'w')
for i in range(np.size(r_nodos)):
	f3.write("{0:22.15e}   {1:22.15e}\n".format(r_nodos[i], ge[i]))
	f4.write("{0:22.15e}   {1:22.15e}\n".format(r_nodos[i], dge[i]))


# exit(0)


for i in range(np.size(r_nodos)):
	f2.write("{0:22.15e}   {1:22.15e}\n".format(a0*r_nodos[i], ge[i]))

f1.close()
f2.close()
