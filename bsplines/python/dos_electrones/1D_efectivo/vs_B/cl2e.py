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
import scipy.sparse.linalg
from scipy import linalg as LA
from bspline import Bspline
from sys import exit
from knots import knots_sequence
from index import indices
from interaccion import V_interaccion
from simetrizacion import H_antisimetrico

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
delta	 = 0.1;
eta		 = 0.1;
B_i		 = 1.0;
B_f 	 = 50.0;

bcampo_vec = np.linspace(B_i, B_f, 49);

# Pendiente de los niveles de landau
slope_ll = eV/(me*c_light*ua_to_T);

## Intervalo de integracion en la coordenada z
Zmax = 1000.0;
Zmin = -Zmax;

N_intervalos_z = 35;
N_cuad = 10;
grado = 6;
kord = grado + 1;
beta = 0.0065; ## Cte de decaimiento para los knots en la distribucion exponencial
N_splines_z = N_intervalos_z + grado;
N_base_z = N_splines_z - 2;

N_dim = N_base_z; # Tamano de las matrices

nev = N_dim-1;

if nev > N_dim:
	print "nev > N_dim, cambiar el valor de nev"
	exit(0)

archivo1 = "./resultados/E_vs_B-v1_%6.4feV-v2_%6.4feV-az_%6.4f-bz_%6.4f-eta_%6.4f.dat" % (v1, v2, az, bz, eta)

f1 = open(archivo1, 'w')
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
f1.write("# Parametro de normalizacion de coulomb 1D delta = {0} en nm\n".format(delta))
f1.write("# Carga efectiva de la interaccion eta = {0} en UA\n".format(eta))
f1.write("# Autovalores calculados\n")

## Paso a unidades atomicas
Zmax = Zmax/a0; Zmin = Zmin/a0;
beta = beta*a0;
az = az/a0; bz = bz/a0; r0 = r0/a0; delta = delta/a0;
v1 = v1/eV; v2 = v2/eV;

## Vector de knots para definir los B-splines, distribucion uniforme
knots_z = knots_sequence(grado, 'exp', N_intervalos_z, beta, Zmin, Zmax);

## Pesos y nodos para la cuadratura de Gauss-Hermite
x, w = np.polynomial.legendre.leggauss(N_cuad);

## Escribo bien los nodos y los pesos para las integrales en z
z_nodos = np.array([]);
w_pesos = np.array([]);
for i in range(N_intervalos_z):
	aux_x = 0.5*(knots_z[i+grado+1]-knots_z[i+grado])*x + 0.5*(knots_z[i+grado+1]+knots_z[i+grado]);
	aux_w = 0.5*(knots_z[i+grado+1]-knots_z[i+grado])*w;

	z_nodos = np.hstack((z_nodos, aux_x));
	w_pesos = np.hstack((w_pesos, aux_w));

w_pesos2 = np.tile(w_pesos, (N_base_z, 1));
z_nodos2 = np.tile(z_nodos, (N_base_z, 1));

## B-splines en la coordenada z
basis = Bspline(knots_z, grado);
## Calculo la matriz con los bsplines y sus derivadas en z
bsz  = [basis._Bspline__basis(i, basis.p) for i in z_nodos]; # evaluo los bsplines
dbsz = [basis.d(i) for i in z_nodos];                        # evaluo las derivadas de los bsplines

bsz = np.array(bsz)
bsz = bsz[:,1:N_splines_z-1]

dbsz = np.array(dbsz)
dbsz = dbsz[:,1:N_splines_z-1]


# Llamo la funcion que calcula la matriz de la interaccion
V_int = V_interaccion(z_nodos, w_pesos, delta, bsz)


### Ahora paso a armar las matrices para el hamiltoniano de una particula
z_mat = np.dot(np.transpose(bsz), (np.transpose(z_nodos2*w_pesos2)*bsz));

## Matriz de solapamiento en z
Sz = np.dot(np.transpose(bsz), (np.transpose(w_pesos2)*bsz));

## Matriz de energia cinetica en z
Tz = 0.5/me*np.dot(np.transpose(dbsz), (np.transpose(w_pesos2)*dbsz));

## Matriz del potencial de confinamiento
Uz = V_potencial_pozo_z1(0.5*az, z_nodos)
Uz = np.tile(Uz, (N_base_z, 1))

Vz1 = np.dot(np.transpose(bsz), (np.transpose(Uz*w_pesos2)*bsz));

Uz = V_potencial_pozo_z2(0.5*az, 0.5*(az+bz), z_nodos)
Uz = np.tile(Uz, (N_base_z, 1))

Vz2 = np.dot(np.transpose(bsz), (np.transpose(Uz*w_pesos2)*bsz));

Tr = 0.5/me*9.74494e-5
Ur = 0.5*me*10677.5
Vr = 0.82604193

## Calculo los autovalores y autovectores
# auval = np.zeros(N_dim);
for bcampo in bcampo_vec:

	l_campo = np.sqrt(alpha/bcampo)/a0;
	omega = 0.5*bcampo/(me*c_light*ua_to_T);  ## Frecuencia de oscilacion debida al campo

	Hr = Tr + omega**2*Ur

	# El hamiltoniano es
	Hz = Tz + v1*Vz2 - v2*Vr*Vz1 + Hr*Sz;

	e1, auvec = LA.eigh(Hz, Sz);
	e1 = eV*e1;

	H_sim, S_sim = H_antisimetrico(Hz, Sz, V_int, eta)

	e, auvec = LA.eig(H_sim, S_sim);
	# e, auvec = scipy.sparse.linalg.eigsh(H_sim, k=np.size(H_sim[0,:])/2, M=S_sim, sigma=None, which='SM')
	e = eV*e;

	# auval = np.vstack((auval, e));

	f1.write("{:13.9e}   ".format(bcampo))
	for i in range(nev):
		f1.write("{:13.9e}   ".format(e.real[i]))

	f1.write("\n")

# for i in range(5):
# 	plt.plot(bcampo_vec, auval[1:,i])
#
# plt.show()
f1.close()
