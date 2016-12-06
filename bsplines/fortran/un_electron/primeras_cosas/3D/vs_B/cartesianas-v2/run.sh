#!/bin/bash

for v in 0.000d0 0.002d0
do
  make V0_=$v

  ./cl1e
  make clean
done
