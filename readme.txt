Problema de confinamiento lateral:

En este repositorio tengo codigos que calculan los autovalores 
para dos electrones confinados por el pozo de potencial:

   V(z) = -V_pozo*exp(-0.5*z^2/sigma^2)

Ademas los electrones interactuan por el potencial de la forma:

  V_int(z1,z2) = sqrt(0.5*pi)/l_pozo*exp(x^2)*(1 - erf(x))

donde x = abs(z1-z2)/(sqrt(2)*l_pozo).
El parametro l_pozo esta dado por la intesidad del campo magnetico, 
l_pozo = sqrt(2*alpha/B), con alpha = 658.4092645439 nm^2/T.

Los programas calculan los autovalores del hamiltoniano de los 
dos electrones usando el metodo variacional de Rayleight-Ritz.

La mejor base que tengo es con la base de los B-splines. Los 
codigos estan escritos en fortran, en C y en octave.
