#include "alhassid_parallel.hpp"   //contains the functions for main program.
#include "extra.hpp"
#include <cstdlib>

int main(int argc, char* argv[])
{
  MPI_Init(NULL, NULL);
  int pRank, num_procs;
  MPI_Comm_size (MPI_COMM_WORLD, &num_procs);
	MPI_Comm_rank (MPI_COMM_WORLD, &pRank);

  if(argc!=5) {cerr << "Enter (1) lattice size, (2) U, (3) fill and (4) no of sweeps.\n"; exit(1);}
  size = atoi(argv[1]);
  U = atof(argv[2]);
  double fill = atof(argv[3]);
  int no_sweeps = atoi(argv[4]);
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

  ofstream outfile_mlength, outfile_data, outfile_dos;
  string dir_path;
  
  if(pRank==0)
  {
    string latticedata = "U_"+to_string(int(U))+"_size_"+to_string(size)+"_sweeps_"+to_string(no_sweeps)+"_t_"+current_time_str();
    dir_path = "pspa_mc/"+latticedata+"/";
    string create_dir_command =  "mkdir -p "+dir_path;
    const int create_dir_result = system(create_dir_command.c_str()); 
	  if(create_dir_result==-1){cerr << "directory creation failed!\n"; exit(1);}

    string filename=dir_path+"data.dat";  outfile_data.open(filename);
    cout << "==============================\n"<< "filename is: " << filename << "\n========================\n";
  }

  vector <double> T_range {0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1,0.08,0.06,0.04,0.02,0.01};
  double final_temp = T_range[0];
  MatrixXcd H_spa = H0-U/2*matrixelement_sigmaz_2d(randsigma);
  pair<MatrixXcd,VectorXd> spa_spectrum = Eigenspectrum(H_spa);
  double spa_F = spa_free_energy(spa_spectrum.second, final_temp, fill);
  double pspa_F = spa_F + pspa_free_energy(final_temp, spa_spectrum.second, spa_spectrum.first.real(), fill);
  
  begin_ms = duration_cast< milliseconds >(system_clock::now().time_since_epoch());

  for(double temperature: T_range)
  {
    for(int sweep=0; sweep< N_therm; sweep++)
    {
      for(int lattice_index=0; lattice_index<size*size; lattice_index++)
      {
        ising_sigma_generate(suggested_randsigma, lattice_index, idum);
        MatrixXcd suggested_Hspa = H0-U/2*matrixelement_sigmaz_2d(suggested_randsigma);
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
      string dos_filename = dir_path+"dos_at_T_"+to_string(temperature)+".dat";
      outfile_dos.open(dos_filename);
      if(outfile_mlength.bad() || outfile_mlength.bad()){cerr << "fail to create file to store aux-fields!"; exit(13);}
    }
    vector <pair <double, double>> dos;
    for(double omega = -U-2*t; omega < U +2*t; omega += 0.1)
    {
      dos.push_back(make_pair(omega,0.0));
    }

    for(int sweep= N_therm; sweep<no_sweeps; sweep++)
    {
      for(int lattice_index=0; lattice_index<size*size; lattice_index++)
      {
        ising_sigma_generate(suggested_randsigma, lattice_index, idum);
        MatrixXcd suggested_Hspa = H0-U/2* matrixelement_sigmaz_2d(suggested_randsigma);
        pair<MatrixXcd,VectorXd> suggested_spa_spectrum = Eigenspectrum(suggested_Hspa);
        double suggested_spa_F = spa_free_energy(suggested_spa_spectrum.second, temperature, fill);
        double suggested_pspa_F = suggested_spa_F 
                                + pspa_free_energy(temperature, suggested_spa_spectrum.second, suggested_spa_spectrum.first.real(), fill);

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

      final_free_energy_spa += spa_F/(size*size);
      final_free_energy_rpa += pspa_F/(size*size); 
      S_pi += get_spi(randsigma);
      if(pRank==0)
      {
        outfile_mlength << sweep <<  " " << randsigma.col(2).transpose() << endl;
        cout << "\r sweep = " << sweep << " done."; cout.flush();
      }
    }
    if(pRank==0)
    {
      for(const auto& it: dos)
      {
        outfile_dos << it.first << " " << it.second/N_meas << endl;
      }
      outfile_data << temperature << " " << final_free_energy_spa/double(N_meas) << " " << final_free_energy_rpa/double(N_meas) << " " << S_pi/double(N_meas) << endl;
      outfile_mlength.close();
      outfile_dos.close();
      cout << "\rtemperature = " << temperature << " done."; cout.flush();
    }
  }

  if(pRank==0)
  {
    cout << endl;
    outfile_data.close();
  }

  MPI_Finalize();
  return 0;
}
