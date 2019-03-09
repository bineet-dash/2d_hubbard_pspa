#include "pspa.hpp"   //contains the functions for main program.
#include "extra.hpp"

int main(int argc, char* argv[])
{
 if(argc!=4) {cerr << "Enter (1) lattice size, (2) U and (3) no of sweeps.\n"; exit(1);}
  size = atoi(argv[1]);
  U = atof(argv[2]);
  int no_sweeps = atoi(argv[3]);
  int N_therm = 0.5*no_sweeps;
  int N_meas = no_sweeps-N_therm;

  int initial_exp = -2;
  int final_exp = -1;
  milliseconds begin_ms, end_ms;
  long idum = time(NULL);

  MatrixXd randsigma=MatrixXd::Zero(size*size,3);
  randsigma.col(2) = VectorXd::Constant(size*size,1);
  for(int i=0; i<randsigma.rows(); i++) ising_sigma_generate(randsigma, i,idum);// randsigma(i,2) = pow(1,xc(i)+yc(i)); //randsigma(i,2) = 5;
  MatrixXd suggested_randsigma = randsigma;
  MatrixXcd H0 = construct_h0_2d(); 
  MatrixXcd Id = MatrixXcd::Identity(H0.rows(),H0.cols());

  ofstream outfile_mlength, outfile_data;
  string dir_path;

  string latticedata = "U_"+to_string(int(U))+"_size_"+to_string(size)+"_sweeps_"+to_string(no_sweeps)+"_t_"+current_time_str();
  dir_path = "spa_mc/"+latticedata+"/";
  string create_dir_command =  "mkdir -p "+dir_path;
  const int create_dir_result = system(create_dir_command.c_str()); 
  if(create_dir_result==-1){cerr << "directory creation failed!\n"; exit(1);}

  string filename=dir_path+"data.dat";  outfile_data.open(filename);
  cout << "==============================\n"<< "filename is: " << filename << "\n========================\n";
  

  vector <double> T_range {0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1,0.08,0.06,0.04,0.02,0.01};
  double final_temp = T_range[0];// 9*pow(10,final_exp);
  MatrixXcd H_spa = H0-U/2*matrixelement_sigmaz_2d(randsigma)+U/4*randsigma.rows()*Id;
  pair<MatrixXcd,VectorXd> spa_spectrum = Eigenspectrum(H_spa);
  double spa_F = spa_free_energy(spa_spectrum.second, final_temp);

  for(double temperature : T_range)
  {
    for(int sweep=0; sweep< N_therm; sweep++)
    {
      for(int lattice_index=0; lattice_index<size*size; lattice_index++)
      {
        ising_sigma_generate(suggested_randsigma, lattice_index, idum);
        MatrixXcd suggested_Hspa = H0-U/2*matrixelement_sigmaz_2d(suggested_randsigma);//+ U/4*randsigma.rows()*Id;
        pair<MatrixXcd,VectorXd> suggested_spa_spectrum = Eigenspectrum(suggested_Hspa);
        double suggested_spa_F = spa_free_energy(suggested_spa_spectrum.second, temperature);

        double uniform_rv = ran0(&idum); 
        double move_prob = exp((spa_F - suggested_spa_F)/temperature);

        if(uniform_rv <= move_prob)
        {
          spa_F = suggested_spa_F;
          randsigma = suggested_randsigma;
        }
        else
        {
          suggested_randsigma=randsigma;
        }
      }
    }

    string fields_filename = dir_path+"field_at_T_"+to_string(temperature)+".dat";
    outfile_mlength.open(fields_filename);
    if(outfile_mlength.bad()){cerr << "fail to create file to store aux-fields!"; exit(13);}
    double final_free_energy = 0; 
    double S_pi =0.0;

    for(int sweep= N_therm; sweep<no_sweeps; sweep++)
    {
      for(int lattice_index=0; lattice_index<size*size; lattice_index++)
      {
        ising_sigma_generate(suggested_randsigma, lattice_index, idum);
        MatrixXcd suggested_Hspa = H0-U/2*matrixelement_sigmaz_2d(suggested_randsigma);//+ U/4*randsigma.rows()*Id;
        pair<MatrixXcd,VectorXd> suggested_spa_spectrum = Eigenspectrum(suggested_Hspa);
        double suggested_spa_F = spa_free_energy(suggested_spa_spectrum.second, temperature);

        double uniform_rv = ran0(&idum); 
        double move_prob = exp((spa_F - suggested_spa_F)/temperature);

        if(uniform_rv <= move_prob)
        {
          spa_F = suggested_spa_F;
          randsigma = suggested_randsigma;
        }
        else
        {
          suggested_randsigma=randsigma;
        }
      }
      final_free_energy += spa_F/(size*size); 
      S_pi += get_spi(randsigma);
      outfile_mlength << sweep <<  " " << randsigma.col(2).transpose() << endl;
      
      cout << "\r sweep = " << sweep << " done."; cout.flush();
    }

    outfile_data << temperature << " " << final_free_energy/N_meas << " " << S_pi/N_meas << " " << endl;
    cout << "\rtemperature = " << temperature << " done."; cout.flush();
    outfile_mlength.close();
  }
 
  cout << endl;
  outfile_data.close();
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

//for(int j=final_exp; j>=initial_exp; j--)
//{
  //for(double i=9; i>=1; i-= decrement)
  //{

  //}