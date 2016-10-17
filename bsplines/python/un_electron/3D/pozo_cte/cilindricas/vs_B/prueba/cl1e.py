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

## Defino el potencial del campo magnetico
def V_campo(me, omega, x):
	U = 0.5*me*omega**2*x**2;
	return U

#############################

## Energia cinetica del momento angular
def L2_angular(me, mz, x):
	L2z = 0.5*mz**2/(me*x**2);
	return L2z

#######################################
#######################################
## Aca empieza el main del programa ###
#######################################
#######################################

## Parametros utiles
a0 = 0.0529177210; eV = 27.21138564; c_light = 137.035999074492; ua_to_T = 1.72e3;

## Parametros fisicos del problema
me		 = 0.063; #0.063; ## Masa de la particula
mz		 = 0.0; ## Componente z del momento angular
B_i		 = 0.0;
B_f 	 = 30.0;

bcampo_vec = np.linspace(B_i, B_f, 30);

## Separo las coordenadas y tomo distinta base en r y en z
Rmin = 0.0;
Rmax = 50.0;

N_intervalos_r = 100;
N_cuad = 200;
grado = 6;
kord = grado + 1;
beta = 0.0065; ## Cte de decaimiento para los knots en la distribucion exponencial
N_splines_r = N_intervalos_r + grado;
N_base_r = N_splines_r - 1;

N_dim = N_base_r; # Tamano de las matrices

# archivo = "./resultados/E_vs_B-v1%6.4feV-v2%6.4feV-Bi%3.1f-Bf%3.1f.dat" % (v1, v2, B_i, B_f)
#
# f = open(archivo, 'w')
# f.write("# Intervalo de integracion en r [{0}, {1}]\n".format(Rmin, Rmax))
# f.write("# Intervalo de integracion en z [{0}, {1}]\n".format(Zmin, Zmax))
# f.write("# Grado de los B-splines {0}\n".format(grado))
# f.write("# Num de intervalos {0} y tamano de base {1} en r\n".format(N_intervalos_r, N_base_r))
# f.write("# Num de intervalos {0} y tamano de base {1} en z\n".format(N_intervalos_z, N_base_z))
# f.write("# Dimension total del espacio N_dim = N_base_r*N_base_z = {0}\n".format(N_dim))
# f.write("# Cte de separacion de knots en dist exp beta = {0}\n".format(beta))
# f.write("# Orden de la cuadratura N_cuad = {0}\n".format(N_cuad))
# f.write("# Masa de la particula me = {0} en UA\n".format(me))
# f.write("# Componente z del momento angular mz = {0}\n".format(mz))
# f.write("# Ancho del pozo en rho rho_0 = {0} en nm\n".format(r0))
# f.write("# Parametros del pozo az = {0} nm y bz = {1} nm\n".format(az, bz))
# f.write("# Altura de la barrera v1 = {0} en eV\n".format(v1))
# f.write("# Profundidad del pozo v2 = {0} en eV\n".format(v2))
# f.write("# Autovalores calculados\n")

## Paso a unidades atomicas
Rmax = Rmax/a0; Rmin = Rmin/a0;
beta = beta*a0;

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
r_nodos2 = np.tile(r_nodos, (N_splines_r, 1));

## B-splines en la coordenada r
basis = Bspline(knots_r, grado);
## Calculo la matriz con los bsplines y sus derivadas en r
bsr  = [basis._Bspline__basis(i, basis.p) for i in r_nodos]; # evaluo los bsplines
dbsr = [basis.d(i) for i in r_nodos];                        # evaluo las derivadas de los bsplines

## Matriz de solapamiento en r
Sr = np.dot(np.transpose(bsr), (np.transpose(r_nodos2)*np.transpose(wr_pesos)*bsr));
Sr = np.array([[Sr[i][j] for i in range(N_splines_r-1)] for j in range(N_splines_r-1)]);

## Matriz de energia cinetica en r
Tr = 0.5/me*np.dot(np.transpose(dbsr), (np.transpose(r_nodos2)*np.transpose(wr_pesos)*dbsr));
Tr = np.array([[Tr[i][j] for i in range(N_splines_r-1)] for j in range(N_splines_r-1)]);

## Matriz de momento angular al cuadrado
L2z = L2_angular(me, mz, r_nodos);
L2z = np.tile(L2z, (N_splines_r, 1));
VL2z = np.dot(np.transpose(bsr), (np.transpose(r_nodos2)*np.transpose(L2z)*np.transpose(wr_pesos)*bsr));
VL2z = np.array([[VL2z[i][j] for i in range(N_splines_r-1)] for j in range(N_splines_r-1)]);

## Matriz del momento angular
Lz = np.dot(np.transpose(bsr), (np.transpose(r_nodos2)*np.transpose(wr_pesos)*bsr));
Lz = np.array([[Lz[i][j] for i in range(N_splines_r-1)] for j in range(N_splines_r-1)]);

## Matriz de energia potencial de campo magnetico
U = V_campo(me, 1.0, r_nodos);
U = np.tile(U, (N_splines_r, 1));
V_B = np.dot(np.transpose(bsr), (np.transpose(r_nodos2)*np.transpose(U)*np.transpose(wr_pesos)*bsr));
V_B = np.array([[V_B[i][j] for i in range(N_splines_r-1)] for j in range(N_splines_r-1)]);

## Calculo los autovalores y autovectores
auval = np.zeros(N_dim);
for bcampo in bcampo_vec:

	omega = 0.5*bcampo/(me*c_light*ua_to_T);  ## Frecuencia de oscilacion debida al campo

	## El hamiltoniano en la coordenada r sin el pozo es
	Hr = Tr + omega**2*V_B + VL2z + omega*mz*Lz;

	e, auvec = LA.eigh(Hr, Sr);
	e = eV*e;

	auval = np.vstack((auval, e));

	# f.write("{:13.9e}   ".format(bcampo))
	# for i in range(10):
	# 	f.write("{:13.9e}   ".format(e[i]))
	#
	# f.write("\n")

for i in range(1):
	plt.plot(bcampo_vec, auval[1:,i])


plt.plot(bcampo_vec, 0.5*1.83758049186099e-3*bcampo_vec)

plt.show()

# f.close()
