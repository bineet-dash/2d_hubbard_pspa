#PBS -N mpiTest
#PBS -l nodes=1:ppn=32
#PBS -j oe 
#PBS -o logs/mpiTest.outlog
#PBS -e logs/mpiTest.errlog
cd /home/bineet/2d_hubbard_pspa 
date
mpiexec -np 2 test.x
date



