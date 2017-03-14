#!/bin/bash

for R in 7.0d0 10.d0 15.d0 20.d0 25.d0 30.d0 35.d0 40.d0 45.d0 50.d0
do
  make RMAX_=$R

  ./cl1e_mod

  make clean
done
