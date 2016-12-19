#!bin/bash/

for carga in 0.001d0 0.0025d0 0.005d0 0.0075d0 0.01d0
do
  make ETA_=$carga

  ./cl2e

  make clean
done
