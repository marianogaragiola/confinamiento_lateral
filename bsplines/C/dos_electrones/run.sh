#!/bin/bash

for z in 50.0 60.0 70.0 80.0 90.0 100.0 110.0 120.0 130.0 140.0 150.0 160.0 170.0 180.0 190.0 200.0 210.0 220.0 230.0
do
  make RMIN=-$z RMAX=$z

  ./cl2e

  make clean
done
