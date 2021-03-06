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
# un polinomio de segundo orden multiplicado por una exponencial
# gaussiana de ancho sigma
#
# V(r,z) = (v1*(r^2+z^2) - v2)*exp(-0.5*(r^2+z^2)/sigma^2)
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
me		 = 0.063; ## Masa de la particula
mz		 = 0.0; ## Componente z del momento angular
sigma	 = 10.0; ## Ancho del pozo gaussiano en nm
v1		 = 0.0042; ## Parametro del pozo en eV/nm^2
v2		 = 0.1; ## Intensidad de campo magnetico en teslas
B_i		 = 0.0;
B_f 	 = 200.0;

bcampo_vec = np.linspace(B_i, B_f, 400);

## Separo las coordenadas y tomo distinta base en r y en z
Rmin = 0.0;
Rmax = 20.0;

Zmax = 50.0;
Zmin = -Zmax;

N_intervalos_r = 30; N_intervalos_z = 30;
N_cuad = 1000;
grado = 6;
kord = grado + 1;
beta = 0.0065; ## Cte de decaimiento para los knots en la distribucion exponencial
N_splines_r = N_intervalos_r + grado; N_splines_z = N_intervalos_z + grado;
N_base_r = N_splines_r - 1; N_base_z = N_splines_z - 2;

N_dim = N_base_r*N_base_z; # Tamano de las matrices

# archivo = "./resultados1909/1e-E_vs_B-v1_%6.4feVsnm2-v2_%6.4feV-sigma_%6.4fnm-mz_%d-Bi_%3.1f-Bf_%3.1f-1.dat" % (v1, v2, sigma, mz, B_i, B_f)
archivo1 = "./resultados1710/E_vs_B-v1%6.4feVsnm2-v2%6.4feV-sigma%6.4fnm-Bi%3.1f-Bf%3.1f-R%6.4f-Z%6.4f.dat" % (v1, v2, sigma, B_i, B_f, Rmax, Zmax)
archivo2 = "./resultados1710/z_vs_B-v1%6.4feVsnm2-v2%6.4feV-sigma%6.4fnm-Bi%3.1f-Bf%3.1f-R%6.4f-Z%6.4f.dat" % (v1, v2, sigma, B_i, B_f, Rmax, Zmax)

f1 = open(archivo1, 'w')
f2 = open(archivo2, 'w')
f1.write("# Intervalo de integracion en r [{0}, {1}]\n".format(Rmin, Rmax))
f1.write("# Intervalo de integracion en z [{0}, {1}]\n".format(Zmin, Zmax))
f1.write("# Grado de los B-splines {0}\n".format(grado))
f1.write("# Num de intervalos {0} y tamano de base {1} en r\n".format(N_intervalos_r, N_base_r))
f1.write("# Num de intervalos {0} y tamano de base {1} en z\n".format(N_intervalos_z, N_base_z))
f1.write("# Dimension total del espacio N_dim = N_base_r*N_base_z = {0}\n".format(N_dim))
f1.write("# Cte de separacion de knots en dist exp beta = {0}\n".format(beta))
f1.write("# Orden de la cuadratura N_cuad = {0}\n".format(N_cuad))
f1.write("# Masa de la particula me = {0} en UA\n".format(me))
f1.write("# Componente z del momento angular mz = {0}\n".format(mz))
f1.write("# Ancho del pozo gaussiano sigma = {0} en nm\n".format(sigma))
f1.write("# Parametro del pozo v1 = {0} en eV/nm^2\n".format(v1))
f1.write("# Parametro del pozo v2 = {0} en eV\n".format(v2))
f1.write("# Autovalores calculados\n")

## Paso a unidades atomicas
Zmax = Zmax/a0; Zmin = Zmin/a0;
Rmax = Rmax/a0; Rmin = Rmin/a0;
sigma = sigma/a0; beta = beta*a0;
v1 = v1*a0**2/eV; v2 = v2/eV;


## Vector de knots para definir los B-splines, distribucion uniforme
knots_r = knots_sequence(grado, 'uniform', N_intervalos_r, beta, Rmin, Rmax);
knots_z = knots_sequence(grado, 'uniform', N_intervalos_z, beta, Zmin, Zmax);

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

### Matrices para el elemento de matriz de z ###
r_mat = np.dot(np.transpose(bsr), (np.transpose(r_nodos2)*np.transpose(wr_pesos)*bsr));
r_mat = np.array([[r_mat[i][j] for i in range(N_splines_r-1)] for j in range(N_splines_r-1)]);

z_mat = np.dot(np.transpose(bsz), (np.transpose(z_nodos2*wz_pesos)*bsz));
z_mat = np.array([[z_mat[i][j] for i in range(1, N_splines_z-1)] for j in range(1, N_splines_z-1)]);

rz_mat = np.kron(r_mat, z_mat)

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
expr = V_potencial_pozo(sigma, r_nodos);
expz = V_potencial_pozo(sigma, z_nodos);

Ur = r_nodos**2*expr;
Uz = z_nodos**2*expz;

expr = np.tile(expr, (N_splines_r, 1));
Ur = np.tile(Ur, (N_splines_r, 1));

expz = np.tile(expz, (N_splines_z, 1));
Uz = np.tile(Uz, (N_splines_z, 1));

V1r = np.dot(np.transpose(bsr), (np.transpose(r_nodos2)*np.transpose(Ur)*np.transpose(wr_pesos)*bsr));
V1r = np.array([[V1r[i][j] for i in range(N_splines_r-1)] for j in range(N_splines_r-1)]);

V2r = np.dot(np.transpose(bsr), (np.transpose(r_nodos2)*np.transpose(expr)*np.transpose(wr_pesos)*bsr));
V2r = np.array([[V2r[i][j] for i in range(N_splines_r-1)] for j in range(N_splines_r-1)]);

V1z = np.dot(np.transpose(bsz), (np.transpose(Uz)*np.transpose(wz_pesos)*bsz));
V1z = np.array([[V1z[i][j] for i in range(1, N_splines_z-1)] for j in range(1, N_splines_z-1)]);

V2z = np.dot(np.transpose(bsz), (np.transpose(expz)*np.transpose(wz_pesos)*bsz));
V2z = np.array([[V2z[i][j] for i in range(1, N_splines_z-1)] for j in range(1, N_splines_z-1)]);

V = v1*(np.kron(V1r, V2z) + np.kron(V2r, V1z)) - v2*np.kron(V2r, V2z);

## Calculo los autovalores y autovectores
auval = np.zeros(N_dim);
for bcampo in bcampo_vec:

	omega = 0.5*bcampo/(me*c_light*ua_to_T);  ## Frecuencia de oscilacion debida al campo

	## El hamiltoniano en la coordenada r sin el pozo es
	Hr = Tr + omega**2*V_B + VL2z + omega*mz*Lz;

	## Armo el hamiltoniano del problema
	Ht = np.kron(Hr, Sz) + np.kron(Sr, Tz) + V;
	St = np.kron(Sr, Sz);

	e, auvec = LA.eigh(Ht, St);
	e = eV*e;

	auval = np.vstack((auval, e));

	z_exp = np.dot(np.transpose(auvec), np.dot(rz_mat, auvec))

	f1.write("{:13.9e}   ".format(bcampo))
	f2.write("{:13.9e}   ".format(bcampo))
	for i in range(30):
		f1.write("{:13.9e}   ".format(e[i]))
		f2.write("{:13.9e}   ".format(z_exp[i,0]))

	f1.write("\n")
	f2.write("\n")


# for i in range(7):
# 	plt.plot(bcampo_vec, auval[1:,i])
#
# plt.show()

f1.close()
f2.close()
