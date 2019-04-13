#!/bin/bash
cd /home/bineet/2d_hubbard_pspa 

#echo "$size"
#echo "$U"
#echo "$nosweeps"

mpiexec -np 8 ./parallel-pspa-ising-mc.x ${size} ${U} ${fill} ${nosweeps}

