#!/bin/bash

for z in 20.0 30.0 40.0 50.0 60.0 70.0 80.0
do
  make B_CAMPO=$z

  ./cl2e

  make clean
done
