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
  milliseconds begin_ms, end_ms;
  long idum = time(NULL);

  MatrixXd randsigma=MatrixXd::Zero(size*size,3);
  for(int i=0; i<randsigma.rows(); i++) randsigma(i,2) = pow(1,xc(i)+yc(i)); //randsigma(i,2) = 5;
  MatrixXd suggested_randsigma = randsigma;
  MatrixXcd H0 = construct_h0_2d(); 
  MatrixXcd Id = MatrixXcd::Identity(H0.rows(),H0.cols());

  MatrixXcd H_spa = H0-U/2*matrixelement_sigmaz_2d(randsigma)+U/4*randsigma.rows()*Id;
  pair<MatrixXcd,VectorXd> spa_spectrum = Eigenspectrum(H_spa);
  pdd free_energies = get_spa_pspa_F(spa_spectrum.first.real(), spa_spectrum.second, final_temp);

  string filename, latticedata;
  latticedata = "U="+to_string(int(U))+"_size="+to_string(size)+"_sweeps="+to_string(no_sweeps);
  filename="pspa_mc/spin_arrangement"+current_time_str()+latticedata+".nb"; ofstream outfile_spinarr(filename);
  spinarrangement_2d_Mathematica_output(randsigma,outfile_spinarr);
  filename="pspa_mc/results_ising_spa_"+current_time_str()+latticedata+".txt"; ofstream outfile_data(filename);
  cout << "==============================\n"<< "filename is: " << filename << "\n========================\n";

  ofstream debugout("debug.txt");
  begin_ms = duration_cast< milliseconds >(system_clock::now().time_since_epoch());

  for(int j=final_exp; j>=initial_exp; j--)
  {
    double decrement = 1;
    // decrement = (j==-1)?0.05:1;
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
          pdd suggested_free_energies = get_spa_pspa_F(suggested_spa_spectrum.first.real(), suggested_spa_spectrum.second, temperature); 
          // debugout << suggested_free_energy << " " << free_energy << endl;

          double uniform_rv = ran0(&idum); double move_prob = exp((free_energies.second - suggested_free_energies.second)/temperature);

          if(uniform_rv <= move_prob)
          {
            free_energies = suggested_free_energies;
            randsigma = suggested_randsigma;
          }
          else
          {
            suggested_randsigma=randsigma;
          }
        }
      }

      double final_free_energy_rpa = 0.0;
      double final_free_energy_spa = 0.0;
      double S_pi = 0.0;

      for(int sweep= N_therm; sweep<no_sweeps; sweep++)
      {
        for(int lattice_index=0; lattice_index<size*size; lattice_index++)
        {
          ising_sigma_generate(suggested_randsigma, lattice_index, idum);
          MatrixXcd suggested_Hspa = H0-U/2* matrixelement_sigmaz_2d(suggested_randsigma)+U/4*randsigma.rows()*Id;
          pair<MatrixXcd,VectorXd> suggested_spa_spectrum = Eigenspectrum(suggested_Hspa);
          pdd suggested_free_energies = get_spa_pspa_F(suggested_spa_spectrum.first.real(), suggested_spa_spectrum.second, temperature); 

          double uniform_rv = ran0(&idum); double move_prob = exp((free_energies.second - suggested_free_energies.second)/temperature);

          if(uniform_rv <= move_prob)
          {
            free_energies = suggested_free_energies;
            randsigma = suggested_randsigma;
          }
          else
          {
            suggested_randsigma=randsigma;
          }
        }
        final_free_energy_spa += free_energies.first;
        final_free_energy_rpa += free_energies.second;

        S_pi += get_spi(randsigma);
        cout << "\r sweep = " << sweep << " done."; cout.flush();
      }

      outfile_data << temperature << " " << final_free_energy_spa/N_meas << " " << final_free_energy_rpa/N_meas << " " << S_pi/N_meas << endl;
      cout << "\rtemperature = " << temperature << " done."; cout.flush();
    }
  }

  cout << endl;
  end_ms = duration_cast< milliseconds >(system_clock::now().time_since_epoch());
  show_time(begin_ms, end_ms,"MC calculation");
  spinarrangement_2d_Mathematica_output(randsigma,outfile_spinarr);

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