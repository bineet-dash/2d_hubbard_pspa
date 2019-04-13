#!/bin/bash
cd /home/bineet/2d_hubbard_pspa 

#echo "$size"
#echo "$U"
#echo "$nosweeps"

mpiexec -np 16 ./parallel-variable-calculation.x ${size} ${U} ${fill} ${inputfile}

