#include "pspa.hpp"
#include "extra.hpp"

milliseconds end_ms, begin_ms;
double Abs(double x){return (x<=0)?-x:x;}

int main(int argc, char* argv[])
{
  if(argc!=5) {cerr << "Enter (1) lattice size, (2) U, (3) fill, (4) dir_path. \n"; exit(1);}
  size = atoi(argv[1]);
  U = atof(argv[2]);
  // double temperature = atof(argv[3]);
  double fill = atof(argv[3]);
  string dir_path = argv[4];

  // cout << "Enter lattice index: ";
  // int li; cin >> li;

  ifstream datain;
  MatrixXd randsigma = MatrixXd::Zero(size*size,3);
  MatrixXcd H0 = construct_h0_2d(); 
  MatrixXcd Id = MatrixXcd::Identity(H0.rows(),H0.cols());

  std::ofstream dataout_s(dir_path+"s"+".dat");
  ofstream dataout_n(dir_path+"n.dat");
  dataout_s.precision(3);
  dataout_n.precision(3);

  vector <double> T_range {0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1,0.08,0.06,0.04,0.02,0.01};
  for(double temperature: T_range)
  {
    VectorXd s_spa = VectorXd::Zero(size*size);
    VectorXd n_spa = VectorXd::Zero(size*size);
    int count = 0;
    string filename = dir_path + "field_at_T_"+to_string(temperature)+".dat";
    datain.open(filename);

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

      vector <MatrixXd> vt_n;
      for(int it=0; it< size*size; it++)
      {
        MatrixXd v_i = MatrixXd::Zero(2*size*size,2*size*size);
        v_i(it,it) = 1; v_i(it+size*size, it+size*size) = 1;
        MatrixXd v_i_transformed = spa_spectrum.first.real().adjoint()*v_i*spa_spectrum.first.real();
        vt_n.push_back(v_i_transformed);
      }

      double mu = get_mu(temperature, spa_spectrum.second, fill);
      VectorXd fermi_hf = VectorXd::Zero(spa_spectrum.second.size());
      for(int it=0; it< spa_spectrum.second.size(); it++)
      {
        fermi_hf(it) = fermi_fn(spa_spectrum.second(it)-mu, temperature);
      }

      for(int li=0; li<size*size; li++)
      {
        for(int i=0; i< fermi_hf.size(); i++) 
        {
          s_spa(li) += vt.at(li)(i,i)*fermi_hf(i); 
          n_spa(li) += vt_n.at(li)(i,i)*fermi_hf(i);
        }
      }
      count ++;
    }
    

    dataout_s << temperature << " " << s_spa.transpose()/double(count) << " " << s_spa.unaryExpr(&Abs).sum()/s_spa.size()/double(count) << endl;
    dataout_n << temperature << " " << n_spa.transpose()/double(count) << endl;
    datain.close();
  }


  
}