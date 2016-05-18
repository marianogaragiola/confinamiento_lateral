#!/bin/bash/

for la in 0.40d0 0.41d0 0.42d0 0.43d0 0.44d0 0.45d0 0.46d0 0.47d0 0.48d0 0.49d0
do 
  make lambda_=$la

  ./2eclvsB

  make clean
done
