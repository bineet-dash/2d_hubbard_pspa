if [ $# -eq 4 ]; then
	qsub -S /bin/bash -N $1 -l nodes=1:ppn=32 -j oe -o logs/$1.outlog -e logs/$1.errlog -v size=$2,U=$3,nosweeps=$4 run_pspa.sh 

else
	echo "Provide 4 arguments: (1) Job name, (2)Size, (3) U, (4) MC sweeps"; exit 10
fi


