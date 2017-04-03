#!bin/bash/

for carga in 0.0d0 0.00001d0 0.0001d0 0.001d0 0.01d0 0.05d0 0.1d0 0.2d0 0.3d0 # 0.4d0 0.5d0 0.6d0 0.7d0 0.8d0 0.9d0 1.0d0
do
  make ETA_=$carga

  ./cl2e

  make clean
done
