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

## Defino el potencial del campo magnetico
def V_campo(me, omega, x):
	U = 0.5*me*omega**2*x**2;
	return U
#############################

## Defino el potencial del pozo
def V_potencial_pozo(sigma, x):
	U = np.exp(-0.5*x**2/sigma**2);
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
me = 0.063; ## Masa de la particula
mz = 0.0; ## Componente z del momento angular
sigma = 15.0; ## Ancho del pozo gaussiano en nm
bcampo = 1.0; ## Intensidad de campo magnetico en teslas

v0_vec = np.linspace(0.0, 0.03, 100);

## Separo las coordenadas y tomo distinta base en r y en z
Rmin = 0.0;
Rmax = 1000.0;

Zmax = 1000.0;
Zmin = -Zmax;

N_intervalos_r = 30; N_intervalos_z = 30;
N_cuad = 100;
grado = 4;
kord = grado + 1;
beta = 0.0065; ## Cte de decaimiento para los knots en la distribucion exponencial
N_splines_r = N_intervalos_r + grado; N_splines_z = N_intervalos_z + grado;
N_base_r = N_splines_r - 1; N_base_z = N_splines_z - 2;

N_dim = N_base_r*N_base_z; # Tamano de las matrices

archivo = "./resultados/1e-E_vs_v0-B_%6.4fT-mz_%d.dat" % (bcampo, mz)

f = open(archivo, 'w')
f.write("# Intervalo de integracion en r [{0}, {1}]\n".format(Rmin, Rmax))
f.write("# Intervalo de integracion en z [{0}, {1}]\n".format(Zmin, Zmax))
f.write("# Grado de los B-splines {0}\n".format(grado))
f.write("# Num de intervalos {0} y tamano de base {1} en r\n".format(N_intervalos_r, N_base_r))
f.write("# Num de intervalos {0} y tamano de base {1} en z\n".format(N_intervalos_z, N_base_z))
f.write("# Dimension total del espacio N_dim = N_base_r*N_base_z = {0}\n".format(N_dim))
f.write("# Cte de separacion de knots en dist exp beta = {0}\n".format(beta))
f.write("# Orden de la cuadratura N_cuad = {0}\n".format(N_cuad))
f.write("# Masa de la particula me = {0} en UA\n".format(me))
f.write("# Componente z del momento angular mz = {0}\n".format(mz))
f.write("# Ancho del pozo gaussiano sigma = {0} en nm\n".format(sigma))
f.write("# Intensidad de campo magnetico B = {0} en T\n".format(bcampo))
f.write("# Autovalores calculados\n")

## Paso a unidades atomicas
Zmax = Zmax/a0; Zmin = Zmin/a0;
Rmax = Rmax/a0; Rmin = Rmin/a0;
sigma = sigma/a0; beta = beta*a0;
omega = 0.5*bcampo/(me*c_light*ua_to_T);  ## Frecuencia de oscilacion debida al campo

## Vector de knots para definir los B-splines, distribucion uniforme
knots_r = knots_sequence(grado, 'exp', N_intervalos_r, beta, Rmin, Rmax);
knots_z = knots_sequence(grado, 'exp', N_intervalos_z, beta, Zmin, Zmax);

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

## Escribo bien los nodos y los pesos para las integrales en z
z_nodos = np.array([]);
wz_pesos = np.array([]);
for i in range(N_intervalos_z+1):
	aux_x = 0.5*(knots_z[i+grado+1]-knots_z[i+grado])*x + 0.5*(knots_z[i+grado+1]+knots_z[i+grado]);
	aux_w = 0.5*(knots_z[i+grado+1]-knots_z[i+grado])*w;

	z_nodos = np.hstack((z_nodos, aux_x));
	wz_pesos = np.hstack((wz_pesos, aux_w));

wz_pesos = np.tile(wz_pesos, (N_splines_z, 1));
z_nodos2 = np.tile(z_nodos, (N_splines_z, 1));


## B-splines en la coordenada r
basis = Bspline(knots_r, grado);
## Calculo la matriz con los bsplines y sus derivadas en r
bsr  = [basis._Bspline__basis(i, basis.p) for i in r_nodos]; # evaluo los bsplines
dbsr = [basis.d(i) for i in r_nodos];                        # evaluo las derivadas de los bsplines

## B-splines en la coordenada z
basis = Bspline(knots_z, grado);
## Calculo la matriz con los bsplines y sus derivadas en z
bsz  = [basis._Bspline__basis(i, basis.p) for i in z_nodos]; # evaluo los bsplines
dbsz = [basis.d(i) for i in z_nodos];                        # evaluo las derivadas de los bsplines


## Matriz de solapamiento en r
Sr = np.dot(np.transpose(bsr), (np.transpose(r_nodos2)*np.transpose(wr_pesos)*bsr));
Sr = np.array([[Sr[i][j] for i in range(N_splines_r-1)] for j in range(N_splines_r-1)]);

## Matriz de solapamiento en z
Sz = np.dot(np.transpose(bsz), (np.transpose(wz_pesos)*bsz));
Sz = np.array([[Sz[i][j] for i in range(1, N_splines_z-1)] for j in range(1, N_splines_z-1)]);

## Matriz de energia cinetica en r
Tr = 0.5/me*np.dot(np.transpose(dbsr), (np.transpose(r_nodos2)*np.transpose(wr_pesos)*dbsr));
Tr = np.array([[Tr[i][j] for i in range(N_splines_r-1)] for j in range(N_splines_r-1)]);

## Matriz de energia cinetica en z
Tz = 0.5/me*np.dot(np.transpose(dbsz), (np.transpose(wz_pesos)*dbsz));
Tz = np.array([[Tz[i][j] for i in range(1, N_splines_z-1)] for j in range(1, N_splines_z-1)]);

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

## Matriz de energia de potencial del pozo de potencial
# Primero en la variable r
U = V_potencial_pozo(sigma, r_nodos);
U = np.tile(U, (N_splines_r, 1));
Vp_r = np.dot(np.transpose(bsr), (np.transpose(r_nodos2)*np.transpose(U)*np.transpose(wr_pesos)*bsr));
Vp_r = np.array([[Vp_r[i][j] for i in range(N_splines_r-1)] for j in range(N_splines_r-1)]);
# Segundo en la variable z
U = V_potencial_pozo(sigma, z_nodos);
U = np.tile(U, (N_splines_z, 1));
Vp_z = np.dot(np.transpose(bsz), (np.transpose(U)*np.transpose(wz_pesos)*bsz));
Vp_z = np.array([[Vp_z[i][j] for i in range(1, N_splines_z-1)] for j in range(1, N_splines_z-1)]);

## El hamiltoniano en la coordenada r es
Hr = Tr + omega**2*V_B + VL2z + omega*Lz;

## Calculo los autovalores y autovectores
auval = np.zeros(N_dim);
for v0 in v0_vec:
	v0 = v0/eV;
	## Armo el hamiltoniano del problema
	Ht = np.kron(Hr, Sz) + np.kron(Sr, Tz) - v0*np.kron(Vp_r, Vp_z);
	St = np.kron(Sr, Sz);

	e, auvec = LA.eigh(Ht, St);
	e = eV*e;

	auval = np.vstack((auval, e));

	f.write("{:13.9e}   ".format(eV*v0))
	for i in range(100):
		f.write("{:13.9e}   ".format(e[i]))

	f.write("\n")

# auval = np.array([[auval[i][j] for i in range(1,np.size(v0_vec)+1)] for j in range(30)]);
#
# for i in range(20):
# 	estado = np.zeros(np.size(v0_vec));
# 	for j in range(np.size(v0_vec)):
# 		estado[j] = auval[i][j];
#
# 	plt.plot(v0_vec, estado, '-')
#
# plt.show()
#
# auval = np.vstack((np.transpose(v0_vec), auval));
#
# np.savetxt(archivo, np.transpose(auval))
f.close()
