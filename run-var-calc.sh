#!/bin/bash
cd /home/bineet/2d_hubbard_pspa 

# echo $size
# echo $U
# echo $T
# echo $fill
# echo $inputfile

mpiexec -np 16 ./debug-parallel-variable-calculation.x ${size} ${U} ${T} ${fill} ${inputfile}

