#!/bin/bash/

for la in 0.3d0 0.4d0 0.5d0 0.6d0 0.7d0
do 
  make lambda_=$la

  ./2eclvsB

  make clean
done
