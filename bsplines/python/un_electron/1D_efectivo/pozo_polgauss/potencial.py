import numpy as np
import matplotlib.pyplot as plt

a0 = 0.0529177210; eV = 27.21138564; c_light = 137.035999074492; ua_to_T = 1.72e3;
alpha = 658.4092645439;

ZMIN = -10.0
ZMAX = 10.0

num_puntos = 100

v1			 = 0.05
v2			 = 0.45
sigma		 = 2.0
b_campo	 = 150.0

archivo = "./potencial/V_QD-v1_%6.4feVsnm2-v2_%6.4feV-B_%6.4fT.dat" % (v1, v2, b_campo)

ZMIN = ZMIN/a0; ZMAX = ZMAX/a0; sigma = sigma/a0
v1 = v1*a0**2/eV; v2 = v2/eV

z = np.linspace(ZMIN, ZMAX, num_puntos)

U = np.exp(-0.5*z**2/sigma**2)

l_campo = np.sqrt(2.0*alpha/b_campo)/a0

v1_2 = v1/(1.0 + l_campo**2/(2.0*sigma**2))
v2_2 = (v1*l_campo**2 - v2)/(1.0 + l_campo**2/(2.0*sigma**2))

V = (v1_2*z**2 + v2_2)*U
V = eV*V

z = a0*z

salida = np.vstack((z, V))

f = open(archivo, 'w')

np.savetxt(archivo, salida)

f.close()

plt.plot(z, V)
plt.show()
