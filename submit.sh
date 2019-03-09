#PBS -N profile_pspa
#PBS -l nodes=1:ppn=4
#PBS -j oe 
#PBS -o logs/out.log
#PBS -e logs/err.log
cd /home/bineet/2d_hubbard_pspa 
date
mpiexec profile-parallel-pspa-ising-mc.x 4 3 0.08  
date
