#!/bin/bash

for R in 5.0d0 6.d0 7.0d0 7.5d0
do
  make R0_=$R AZ_=$R RMAX_=$R

  ./cl1e

  make clean
done
