#!/bin/bash
#
#$ -cwd
#
#$ -q blcr@compute-0-13.local
#
#$ -m eas -M marianogaragiola@gmail.com

make run

make clean
