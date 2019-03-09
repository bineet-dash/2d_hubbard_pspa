#include "alhassid_parallel.hpp"   //contains the functions for main program.
#include "extra.hpp"
#include <cstdlib>

int main(int argc, char* argv[])
{
  MPI_Init(NULL, NULL);
  int pRank, num_procs;
  MPI_Comm_size (MPI_COMM_WORLD, &num_procs);
	MPI_Comm_rank (MPI_COMM_WORLD, &pRank);

  if(argc!=4) {cerr << "Enter (1) lattice size, (2) U and (3) temperature.\n"; exit(1);}
  size = atoi(argv[1]);
  U = atof(argv[2]);
  double temperature = atof(argv[3]);

  long idum = time(NULL);
  milliseconds begin_ms, end_ms;
  begin_ms = duration_cast< milliseconds >(system_clock::now().time_since_epoch());

  MatrixXd randsigma=MatrixXd::Zero(size*size,3);
  randsigma.col(2) = VectorXd::Constant(size*size,1);
  for(int i=0; i<randsigma.rows(); i++) randsigma(i,2) = pow(1,xc(i)+yc(i));
  MatrixXd suggested_randsigma = randsigma;
  MatrixXcd H0 = construct_h0_2d(); 
  MatrixXcd Id = MatrixXcd::Identity(H0.rows(),H0.cols());

  MatrixXcd H_spa = H0-U/2*matrixelement_sigmaz_2d(randsigma);
  pair<MatrixXcd,VectorXd> spa_spectrum = Eigenspectrum(H_spa);
  
  double spa_F = spa_free_energy(spa_spectrum.second, temperature);
  double pspa_F = spa_F + profile_pspa_free_energy(temperature, spa_spectrum.second, spa_spectrum.first.real(), pRank);
  end_ms = duration_cast< milliseconds >(system_clock::now().time_since_epoch());
  
  if(pRank==0)
  {
    cout << spa_F/(size*size) << " " << pspa_F/(size*size) << endl;
    cout << "Code took " <<  double((end_ms-begin_ms).count())/1000.0 << " seconds. " << endl; 
  }

  MPI_Finalize();
  return 0;
}

/* vector <double> debug;
for(int ix=0; ix<size; ix++)
{
  for(int iy=0; iy<size; iy++)
  {
    double kx = 2*M_PI*ix/size-M_PI;
    double ky = 2*M_PI*iy/size-M_PI;
    debug.push_back(filter(-2*t*(cos(kx)+cos(ky))));
  }
}
sort(debug.begin(),debug.end());  
for(auto it=debug.begin(); it!=debug.end(); it++) cout << filter(*it) << " "; cout << endl;
cout << endl; exit(1); */