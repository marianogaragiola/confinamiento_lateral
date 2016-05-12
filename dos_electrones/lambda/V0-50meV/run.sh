
#!/bin/bash
#
#$ -cwd
#
for B in 5 10 20 40 60
do

make B_campo_=$B

./2eclvslambda

make clean
done
