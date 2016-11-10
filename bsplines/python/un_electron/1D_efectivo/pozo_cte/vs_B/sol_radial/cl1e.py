# Codigo que calcula los autovalores y autovectores del hamiltoniano
# efectio 1D para un electron.
#
# El hamiltoniano del sistema es
#
#   Hz = -hbar^2/(2*me)*d^2/dz^2 + V_{ef}(z)
#
# El potencial de confinamiento efectivo es
#
#   V_{ef}(z) = (v1*l^2/(1+l^2/(2*sigma^2))^2 + (v1*z^2-v2)/(1+l^2/(2*sigma^2)))*exp(-z^2/(2*sigma^2))
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

## Defino el potencial del pozo en la coordenada z
def V_potencial_pozo_z1(z1, x):
	U = np.piecewise(x, [x<-z1, (-z1<=x)&(x<z1), z1<x], [0, 1, 0]);
	return U

def V_potencial_pozo_z2(z1, z2, x):
	U = np.piecewise(x, [x<-z2, (-z2<=x)&(x<-z1), (-z1<=x)&(x<=z1), (z1<x)&(x<=z2), z2<x], [0, 1, 0, 1, 0]);
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
me		 = 0.041; ## Masa de la particula
r0		 = 7.0; ## Ancho del pozo gaussiano en nm
az		 = 7.0;
bz		 = 2.5;
v1		 = 0.37; ## Alto de la barrera
v2		 = 0.108844; ## Profundidad del pozo
B_i		 = 1.0;
B_f 	 = 50.0;

bcampo_vec = np.linspace(B_i, B_f, 60);

# Pendiente de los niveles de landau
slope_ll = eV/(me*c_light*ua_to_T);

## Intervalo de integracion en la coordenada z
Zmax = 1000.0;
Zmin = -Zmax;

N_intervalos_z = 100;
N_cuad = 100;
grado = 6;
kord = grado + 1;
beta = 0.0065; ## Cte de decaimiento para los knots en la distribucion exponencial
N_splines_z = N_intervalos_z + grado;
N_base_z = N_splines_z - 2;

N_dim = N_base_z; # Tamano de las matrices

nev = N_dim;

if nev > N_dim:
	print "nev > N_dim, cambiar el valor de nev"
	exit(0)

archivo1 = "./resultados/E_vs_B-v1_%6.4feV-v2_%6.4feV-az_%6.4f-bz_%6.4f.dat" % (v1, v2, az, bz)
archivo2 = "./resultados/z_vs_B-v1_%6.4feV-v2_%6.4feV-az_%6.4f-bz_%6.4f.dat" % (v1, v2, az, bz)

f1 = open(archivo1, 'w')
f2 = open(archivo2, 'w')
f1.write("# Intervalo de integracion en z [{0}, {1}]\n".format(Zmin, Zmax))
f1.write("# Grado de los B-splines {0}\n".format(grado))
f1.write("# Num de intervalos {0} y tamano de base {1} en z\n".format(N_intervalos_z, N_base_z))
f1.write("# Dimension total del espacio N_dim = N_base_z = {0}\n".format(N_dim))
f1.write("# Cte de separacion de knots en dist exp beta = {0}\n".format(beta))
f1.write("# Orden de la cuadratura N_cuad = {0}\n".format(N_cuad))
f1.write("# Masa de la particula me = {0} en UA\n".format(me))
f1.write("# Parametro del pozo az = {0} en nm\n".format(az))
f1.write("# Parametro del pozo bz = {0} en nm\n".format(bz))
f1.write("# Parametro del pozo r0 = {0} en nm\n".format(r0))
f1.write("# Parametro del pozo v1 = {0} en eV\n".format(v1))
f1.write("# Parametro del pozo v2 = {0} en eV\n".format(v2))
f1.write("# Autovalores calculados\n")

## Paso a unidades atomicas
Zmax = Zmax/a0; Zmin = Zmin/a0;
beta = beta*a0;
az = az/a0; bz = bz/a0; r0 = r0/a0;
v1 = v1/eV; v2 = v2/eV;

## Vector de knots para definir los B-splines, distribucion uniforme
knots_z = knots_sequence(grado, 'exp', N_intervalos_z, beta, Zmin, Zmax);

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
z_nodos2 = np.tile(z_nodos, (N_splines_z, 1));

## B-splines en la coordenada z
basis = Bspline(knots_z, grado);
## Calculo la matriz con los bsplines y sus derivadas en z
bsz  = [basis._Bspline__basis(i, basis.p) for i in z_nodos]; # evaluo los bsplines
dbsz = [basis.d(i) for i in z_nodos];                        # evaluo las derivadas de los bsplines

z_mat = np.dot(np.transpose(bsz), (np.transpose(z_nodos2*wz_pesos)*bsz));
z_mat = np.array([[z_mat[i][j] for i in range(1, N_splines_z-1)] for j in range(1, N_splines_z-1)]);

## Matriz de solapamiento en z
Sz = np.dot(np.transpose(bsz), (np.transpose(wz_pesos)*bsz));
Sz = np.array([[Sz[i][j] for i in range(1, N_splines_z-1)] for j in range(1, N_splines_z-1)]);

## Matriz de energia cinetica en z
Tz = 0.5/me*np.dot(np.transpose(dbsz), (np.transpose(wz_pesos)*dbsz));
Tz = np.array([[Tz[i][j] for i in range(1, N_splines_z-1)] for j in range(1, N_splines_z-1)]);

## Matriz del potencial de confinamiento
Uz = V_potencial_pozo_z1(0.5*az, z_nodos)
Uz = np.tile(Uz, (N_splines_z, 1))

Vz1 = np.dot(np.transpose(bsz), (np.transpose(Uz*wz_pesos)*bsz));
Vz1 = np.array([[Vz1[i][j] for i in range(1, N_splines_z-1)] for j in range(1, N_splines_z-1)]);

Uz = V_potencial_pozo_z2(0.5*az, 0.5*(az+bz), z_nodos)
Uz = np.tile(Uz, (N_splines_z, 1))

Vz2 = np.dot(np.transpose(bsz), (np.transpose(Uz*wz_pesos)*bsz));
Vz2 = np.array([[Vz2[i][j] for i in range(1, N_splines_z-1)] for j in range(1, N_splines_z-1)]);

Tr = 0.5/me*9.74494e-5
Ur = 0.5*me*10677.5
Vr = 0.82604193

## Calculo los autovalores y autovectores
auval = np.zeros(N_dim);
for bcampo in bcampo_vec:

	l_campo = np.sqrt(alpha/bcampo)/a0;
	omega = 0.5*bcampo/(me*c_light*ua_to_T);  ## Frecuencia de oscilacion debida al campo

	Hr = Tr + omega**2*Ur

	# El hamiltoniano es

	Hz = Tz + v1*Vz2 - v2*Vr*Vz1 + Hr*Sz;

	e, auvec = LA.eigh(Hz, Sz);
	e = eV*e;

	# e = e + 0.5*slope_ll*bcampo

	auval = np.vstack((auval, e));

	z_exp = np.dot(np.transpose(auvec), np.dot(z_mat, auvec))

	f1.write("{:13.9e}   ".format(bcampo))
	f2.write("{:13.9e}   ".format(bcampo))
	for i in range(nev):
		f1.write("{:13.9e}   ".format(e[i]))
		f2.write("{:13.9e}   ".format(z_exp[i,0]))

	f1.write("\n")
	f2.write("\n")

for i in range(20):
	plt.plot(bcampo_vec, auval[1:,i])

plt.show()
f1.close()
f2.close()
