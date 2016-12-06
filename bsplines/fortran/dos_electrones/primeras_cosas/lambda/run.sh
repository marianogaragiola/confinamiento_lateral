#!/bin/bash

for z in 50.d0 60.d0 70.d0 80.d0 90.d0 100.d0 110.d0 120.d0 130.d0 140.d0 150.d0 160.d0 170.d0 180.d0 190.d0 200.d0 210.d0 220.d0 230.d0
do
  make zmin_=-$z zmax_=$z

  ./2eclvslambda

  make clean
done
