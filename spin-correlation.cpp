#include "pspa.hpp"
#include "extra.hpp"

milliseconds end_ms, begin_ms;

double corr_ri(int i, int r, const MatrixXd& spin_lattice)
{
  double corr = 0;
  for(int yd = 1; yd < r; yd++)
  {
    int xd = r-yd;
    corr += spin_lattice(xc(i),yc(i))*spin_lattice(periodic(xc(i),xd, size), periodic(yc(i), yd, size));
  }
  for(int xd = 1; xd < r; xd++)
  {
    int yd = r-xd;
    corr += spin_lattice(xc(i),yc(i))*spin_lattice(periodic(xc(i),xd, size), periodic(yc(i), yd, size));
  }
  return corr;
}

double corr_r(int r, const MatrixXd& spin_lattice)
{
  double corr = 0;
  for(int i=0; i<size*size; i++)
  {
    corr += corr_ri(i, r, spin_lattice);
  }
  return corr/(size*size);
}

int main(int argc, char* argv[])
{
  if(argc!=3) {cerr << "Enter (1) lattice size, (3) spin dir_path. \n"; exit(1);}
  size = atoi(argv[1]);
  string dir_path = argv[2];
  string file_path = dir_path+"s.dat";
  ifstream config_in (file_path);

  ofstream dataout;
  // vector <double> T_range {0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1,0.08,0.06,0.04,0.02,0.01};
  while(true)
  {
    double temperature; 
    config_in >> temperature;
    string filename = dir_path + "spin_corr_at"+to_string(temperature)+".dat";
    dataout.open(filename);
    
    MatrixXd spin_lattice = MatrixXd::Zero(size,size);
    for(int i=0; i<size*size; i++)
    {
      config_in >> spin_lattice(xc(i), yc(i));
    }
    double avg_moment;
    config_in >> avg_moment;
    
    for(int r=0; r<size; r++)
    {
      dataout << r << " " << corr_r(r, spin_lattice) << endl;
    }
    dataout.close();

    if(temperature==0.5)
    {
      cout << "spins: " << endl;
      cout << spin_lattice << endl << endl;
    } 
    if(config_in.fail()) break;
  }

}