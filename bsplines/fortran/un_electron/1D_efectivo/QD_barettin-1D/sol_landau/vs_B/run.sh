#!/bin/bash

for R in 1.d0 2.d0 3.d0 4.d0 5.d0
do
  make R0_=$R RMAX_=$R

  ./cl1e

  make clean
done
