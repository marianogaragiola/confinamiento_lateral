#!bin/bash/

for carga in 0.0075d0 0.01d0
do
  make ETA_=$carga

  ./cl2e

  make clean
done
