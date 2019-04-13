#include "pspa.hpp"
#include "extra.hpp"


milliseconds end_ms, begin_ms;

int main(int argc, char* argv[])
{
  if(argc!=6) {cerr << "Enter (1) lattice size, (2) U, (3) temperature, (4) fill, (5) input filename. \n"; exit(1);}
  size = atoi(argv[1]);
  U = atof(argv[2]);
  double temperature = atof(argv[3]);
  double fill = atof(argv[4]);

  ifstream datain;
  datain.open(argv[4]);

  MatrixXd randsigma = MatrixXd::Zero(size*size,3);
  MatrixXcd H0 = construct_h0_2d(); 
  MatrixXcd Id = MatrixXcd::Identity(H0.rows(),H0.cols());


  double d_spa = 0;
  int count = 0;
  while(!datain.eof())
  {
    int sweep; datain >> sweep;
    for(int i=0; i<size*size; i++) datain >> randsigma(i, 2);
    MatrixXcd H_spa = H0-U/2*matrixelement_sigmaz_2d(randsigma);//+U/4*randsigma.rows()*Id;
    pair<MatrixXcd,VectorXd> spa_spectrum = Eigenspectrum(H_spa);
    
    vector <MatrixXd> vt;
    for(int it=0; it< size*size; it++)
    {
      MatrixXd v_i = MatrixXd::Zero(2*size*size,2*size*size);
      v_i(it,it) = 1; v_i(it+size*size, it+size*size) = -1;
      MatrixXd v_i_transformed = spa_spectrum.first.real().adjoint()*v_i*spa_spectrum.first.real();
      vt.push_back(v_i_transformed);
    }

    double mu = get_mu(temperature, spa_spectrum.second, fill);
    VectorXd fermi_hf = VectorXd::Zero(spa_spectrum.second.size());
    for(int it=0; it< spa_spectrum.second.size(); it++)
    {
      fermi_hf(it) = fermi_fn(spa_spectrum.second(it)-mu, temperature);
    }

    for(int i=0; fermi_hf.size(); i++) d_spa += vt.at(0)(i,i)*fermi_hf(i); 
    count ++;
  }

  cout << d_spa << endl;
  
}