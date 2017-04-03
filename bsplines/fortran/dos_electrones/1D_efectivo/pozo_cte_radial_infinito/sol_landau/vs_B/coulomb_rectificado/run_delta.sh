#!bin/bash/

for d in 0.0001d0 0.0005d0 0.001d0 0.005d0 0.01d0 0.05d0 0.1d0
do
  make DELTA_=$d

  ./cl2e

  make clean
done
