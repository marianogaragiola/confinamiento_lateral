#!/bin/bash

for i in 0.d0  #0.1d0 0.2d0 0.3d0 0.4d0 0.5d0 0.6d0 0.7d0 0.9d0 0.9d0 1.0d0
do
  gfortran -cpp -O3 -o cl2e4 cl2e4.f90 -llapack -Dme_=0.063d0 -DV0_=0.05d0 -DB_campo_=80.d0 -Dlambda_=$i -Domega_i_=0.000001d0 -Domega_f_=0.03d0 -DN_omega_=100 -Dsigma_=20.d0

  ./cl2e4

  rm cl2e4 *.mod *.o
done
