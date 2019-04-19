if [ $# -eq 6 ]; then
	qsub -S /bin/bash -N $1 -l nodes=1:ppn=32 -j oe -o logs/$1.outlog -e logs/$1.errlog -v size=$2,U=$3,fill=$4,T=$5,inputfile=$6 run-var-calc.sh 

else
	echo "Provide 6 arguments: (1) Job name, (2)Size, (3) U, (4) fill, (5) temperature, (6) inputfile"; exit 10
fi


