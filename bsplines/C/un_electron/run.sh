#!/bin/bash/

for b in 0.004 0.0045 0.005 0.0055 0.006 0.0065 0.007 0.008 0.009 0.01
do 
	make BETA=$b

	./cl1e

	make clean
done
