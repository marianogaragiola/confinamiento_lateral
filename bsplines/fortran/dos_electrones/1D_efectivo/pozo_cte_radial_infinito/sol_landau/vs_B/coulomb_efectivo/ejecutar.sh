#!bin/bash/

for carga in 0.0d0 0.00001d0 0.00002d0 0.00005d0 0.00007d0 0.0001d0 0.001d0 0.01d0
do
  make ETA_=$carga

  ./cl2e

  make clean
done
