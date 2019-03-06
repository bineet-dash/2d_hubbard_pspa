#include "alhassid_parallel.hpp"   //contains the functions for main program.
#include "extra.hpp"
#include <cstdlib>

int main(int argc, char* argv[])
{
  MPI_Init(NULL, NULL);
  int pRank, num_procs;
  MPI_Comm_size (MPI_COMM_WORLD, &num_procs);
	MPI_Comm_rank (MPI_COMM_WORLD, &pRank);

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
  
  if(pRank==0)
  {
    string latticedata = "_U_"+to_string(int(U))+"_size_"+to_string(size)+"_sweeps_"+to_string(no_sweeps)+"_t_"+current_time_str();
    dir_path = "pspa_mc/"+latticedata+"/";
    string create_dir_command =  "mkdir -p "+dir_path;
    const int create_dir_result = system(create_dir_command.c_str()); 
	  if(create_dir_result==-1){cerr << "directory creation failed!\n"; exit(1);}

    string filename=dir_path+"data.dat";  outfile_data.open(filename);
    cout << "==============================\n"<< "filename is: " << filename << "\n========================\n";
  }

  double final_temp = 9*pow(10,final_exp);
  MatrixXcd H_spa = H0-U/2*matrixelement_sigmaz_2d(randsigma)+U/4*randsigma.rows()*Id;
  pair<MatrixXcd,VectorXd> spa_spectrum = Eigenspectrum(H_spa);
  double spa_F = spa_free_energy(spa_spectrum.second, final_temp);
  double pspa_F = spa_F + pspa_free_energy(final_temp, spa_spectrum.second, spa_spectrum.first.real());
  
  begin_ms = duration_cast< milliseconds >(system_clock::now().time_since_epoch());

  for(int j=final_exp; j>=initial_exp; j--)
  {
    double decrement = 1;
    for(double i=9; i>=1; i-= decrement)
    {
      double temperature = i*pow(10,j);
      for(int sweep=0; sweep< N_therm; sweep++)
      {
        for(int lattice_index=0; lattice_index<size*size; lattice_index++)
        {
          ising_sigma_generate(suggested_randsigma, lattice_index, idum);
          MatrixXcd suggested_Hspa = H0-U/2*matrixelement_sigmaz_2d(suggested_randsigma)+ U/4*randsigma.rows()*Id;
          pair<MatrixXcd,VectorXd> suggested_spa_spectrum = Eigenspectrum(suggested_Hspa);
          double suggested_spa_F = spa_free_energy(suggested_spa_spectrum.second, temperature);
          double suggested_pspa_F = suggested_spa_F + pspa_free_energy(temperature, suggested_spa_spectrum.second, suggested_spa_spectrum.first.real());

          double move_prob = exp(-(suggested_pspa_F-pspa_F)/temperature);
          double uniform_rv = ran0(&idum);

          if(uniform_rv <= move_prob)
          {
            randsigma = suggested_randsigma;
            spa_F = suggested_spa_F;
            pspa_F = suggested_pspa_F;
          }
          else
          {
            suggested_randsigma = randsigma;
          }
        }
        if(pRank==0) cout << "\r sweep = " << sweep << " done."; cout.flush();
      }

      double final_free_energy_rpa = 0.0;
      double final_free_energy_spa = 0.0;
      double S_pi = 0.0;

      if(pRank==0)
      {
        string filename = dir_path+"field_at_T_"+to_string(temperature)+".dat";
        outfile_mlength.open(filename);
      }

      for(int sweep= N_therm; sweep<no_sweeps; sweep++)
      {
        for(int lattice_index=0; lattice_index<size*size; lattice_index++)
        {
          ising_sigma_generate(suggested_randsigma, lattice_index, idum);
          MatrixXcd suggested_Hspa = H0-U/2* matrixelement_sigmaz_2d(suggested_randsigma)+U/4*randsigma.rows()*Id;
          pair<MatrixXcd,VectorXd> suggested_spa_spectrum = Eigenspectrum(suggested_Hspa);
          double suggested_spa_F = spa_free_energy(suggested_spa_spectrum.second, temperature);
          double suggested_pspa_F = suggested_spa_F + pspa_free_energy(temperature, suggested_spa_spectrum.second, suggested_spa_spectrum.first.real());

          double move_prob = exp(-(suggested_pspa_F-pspa_F)/temperature);
          double uniform_rv = ran0(&idum);

          if(uniform_rv <= move_prob)
          {
            randsigma = suggested_randsigma;
            spa_F = suggested_spa_F;
            pspa_F = suggested_pspa_F;
          }
          else
          {
            suggested_randsigma = randsigma;
          }

          final_free_energy_spa += spa_F/(size*size);
          final_free_energy_rpa += pspa_F/(size*size); 
          S_pi += get_spi(randsigma);
          if(pRank==0)
          {
            outfile_mlength << temperature <<  " " << randsigma.col(2).transpose() << endl;
            cout << "\r sweep = " << sweep << " done."; cout.flush();
          }
        }
      }
      if(pRank==0)
      {
        outfile_data << temperature << " " << final_free_energy_spa/double(N_meas) << " " << final_free_energy_rpa/double(N_meas) << " " << S_pi/double(N_meas) << endl;
        cout << "\rtemperature = " << temperature << " done."; cout.flush();
      }
    }
  }

  if(pRank==0)
  {
    cout << endl;
    outfile_mlength.close();
    outfile_data.close();
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