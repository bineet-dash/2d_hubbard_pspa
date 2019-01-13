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

  int initial_exp = -3;
  int final_exp = 0;
  double final_temp = 10*pow(10,final_exp);

  milliseconds begin_ms = std::chrono::duration_cast< milliseconds >(system_clock::now().time_since_epoch());

  MatrixXd randsigma=MatrixXd::Zero(size*size,3);
  long idum = time(NULL);
  for(int i=0; i<randsigma.rows(); i++) randsigma(i,2) = pow(1,xc(i)+yc(i)); //randsigma(i,2) = 5;

  MatrixXcd H0 = construct_h0_2d(); 
  MatrixXcd Id = MatrixXcd::Identity(H0.rows(),H0.cols());

  MatrixXcd initial_Hamiltonian = H0-U/2*matrixelement_sigmaz_2d(randsigma)+U/4*randsigma.rows()*Id;
  // cout << initial_Hamiltonian << endl << endl << endl;
  double internal_energy = spa_internal_energy(initial_Hamiltonian, final_temp);

  MatrixXd suggested_randsigma = randsigma;

  string filename, latticedata;
  latticedata = "_U_"+to_string(int(U))+"_size_"+to_string(size)+"_sweeps_"+to_string(no_sweeps);
  filename="spa_mc/spin_arrangement_internal_energy_"+current_time_str()+latticedata+".nb"; ofstream outfile_spinarr(filename);
  spinarrangement_2d_Mathematica_output(randsigma,outfile_spinarr);
  filename="spa_mc/results_ising_spa_internal_energy_"+current_time_str()+latticedata+".txt"; ofstream outfile_data(filename);
  cout << "==============================\n"<< "filename is: " << filename << "\n========================\n";

  // for(int j=final_exp; j>=initial_exp; j--)
  // {
  //   double decrement = 0.1;
  //   // decrement = (j==-1)?0.05:1;
  //   for(double i=9; i>=1; i-= decrement)
  //   {
  //     double temperature = i*pow(10,j);
    for(double temperature = 2.0; temperature >= 0.01; temperature -= 0.01)
    {
      for(int sweep=0; sweep< N_therm; sweep++)
      {
        for(int lattice_index=0; lattice_index<size*size; lattice_index++)
        {
          ising_sigma_generate(suggested_randsigma, lattice_index, idum);
          MatrixXcd suggested_Hamiltonian = H0-U/2*matrixelement_sigmaz_2d(suggested_randsigma)+ U/4*randsigma.rows()*Id;
          double suggested_internal_energy = spa_internal_energy(suggested_Hamiltonian, temperature);

          double uniform_rv = ran0(&idum); double move_prob = exp((internal_energy - suggested_internal_energy)/temperature);

          if(uniform_rv <= move_prob)
          {
            internal_energy = suggested_internal_energy;
            randsigma = suggested_randsigma;
          }
          else
          {
            suggested_randsigma=randsigma;
          }
        }
      }

      double final_internal_energy = 0; double S_pi =0.0; double magnetisation=0.0; 
      for(int sweep= N_therm; sweep<no_sweeps; sweep++)
      {
        for(int lattice_index=0; lattice_index<size*size; lattice_index++)
        {
          ising_sigma_generate(suggested_randsigma, lattice_index, idum);
          MatrixXcd suggested_Hamiltonian = H0-U/2* matrixelement_sigmaz_2d(suggested_randsigma)+U/4*randsigma.rows()*Id;
          double suggested_internal_energy = spa_internal_energy(suggested_Hamiltonian, temperature);

          double uniform_rv = ran0(&idum); double move_prob = exp((internal_energy - suggested_internal_energy)/temperature);

          if(uniform_rv <= move_prob)
          {
            internal_energy = suggested_internal_energy;
            randsigma = suggested_randsigma;
          }
          else
          {
            suggested_randsigma=randsigma;
          }
        }
        final_internal_energy += internal_energy/size*size; 

        S_pi += get_spi(randsigma);
        cout << "\r sweep = " << sweep << " done."; cout.flush();
      }

      outfile_data << temperature << " " << final_internal_energy/N_meas << " " << S_pi/N_meas << endl;
      cout << "\rtemperature = " << temperature << " done."; cout.flush();
    }
  // }

  cout << endl;
  milliseconds end_ms = duration_cast< milliseconds >(system_clock::now().time_since_epoch());
  show_time(begin_ms, end_ms,"MC calculation");
  
  // outfile_mcdetails.close();
  outfile_data.close();
  outfile_spinarr.close();
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