#include "pspa.hpp"
#include "extra.hpp"
#include "alhassid_parallel.hpp"


int main(int argc, char* argv[])
{
  if(argc!=5) {cerr << "Enter (1) lattice size, (2) U, (3) temperature, (4) input filename. \n"; exit(1);}
  size = atoi(argv[1]);
  U = atof(argv[2]);
  double temperature = atof(argv[3]);

  ifstream datain;
  datain.open(argv[4]);

  MatrixXd randsigma = MatrixXd::Zero(size*size,3);
  MatrixXcd H0 = construct_h0_2d(); 
  MatrixXcd Id = MatrixXcd::Identity(H0.rows(),H0.cols());


  // while(!datain.eof())
  // {
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



    double mu = get_mu(temperature, spa_spectrum.second);
    VectorXd fermi_hf = VectorXd::Zero(spa_spectrum.second.size());
    for(int it=0; it< spa_spectrum.second.size(); it++)
    {
      fermi_hf(it) = fermi_fn(spa_spectrum.second(it)-mu, temperature);
    }

    cd d_pspa = 0.0;

    // cout << spa_spectrum.second.size() << endl;

    for(int Is = 0; Is < size*size; Is++)
    {
      for(int Js = 0; Js < size*size; Js++)
      {
        
        cd*** I_ijk = new cd** [2*size];
        for(int i=0; i<2*size; i++) I_ijk[i] = new cd* [2*size];
        for(int i=0; i<2*size; i++)
        {
          for(int j=0; j<2*size; j++) I_ijk[i][j] = new cd [2*size];
        }

        // for(int i=0; i<2*size; i++)
        // {
        //   for(int j=0; j<2*size; j++)
        //   {
        //     for(int k=0; k<2*size; k++)
        //     {
        //       I_ijk[i][j][k] = 0.0; 
        //     }
        //   } 
        // }

        milliseconds begin_ms = duration_cast< milliseconds >(system_clock::now().time_since_epoch());
        for(int i=0; i<2*size; i++)
        {
          for(int j=0; j<2*size; j++)
          {
            for(int k=0; k<2*size; k++)
            {
              for(int r=0; r<100; r++)
              {
                cd temp = 0;
                temp += I_ijkr(i,j,k,r,temperature, fermi_hf, spa_spectrum.second)
                                  /(1.0 + a_IJr(Is,Js,r,temperature, spa_spectrum.second, fermi_hf, vt));
                I_ijk[i][j][k] = temp;
              }
            }
          } 
        }
        milliseconds end_ms = duration_cast < milliseconds> (system_clock::now().time_since_epoch());
        cout << end_ms.count()- begin_ms.count() << endl;

        begin_ms = duration_cast< milliseconds >(system_clock::now().time_since_epoch());
        MatrixXcd B_ki = MatrixXcd::Zero(spa_spectrum.second.size(), spa_spectrum.second.size());

        for(int i=0; i<spa_spectrum.second.size(); i++)
        {
          // cout << "i = " << i << " started. \n";
          for(int k=0; k<spa_spectrum.second.size(); k++)
          {
            for(int j=0; j<spa_spectrum.second.size(); j++)
            {
              B_ki(i,k) += vt.at(Is)(i,j)*vt.at(Js)(j,k)*vt.at(0)(k,i)*I_ijk[i][j][k];
            }
          }
        }
        delete[] I_ijk;

        end_ms = duration_cast < milliseconds> (system_clock::now().time_since_epoch());
        cout << end_ms.count()- begin_ms.count() << endl;
        exit(1);

      }
    }

    cout << d_pspa << endl;
    // ofstream varout("ar.dat");
    // for(int r=0; r<1000; r++)
    // {
    //   varout << abs(a_IJr(0 ,1, r, temperature, spa_spectrum.second, fermi_hf, vt)) << " \t " << abs(I_ijkr(0,2,4,r,temperature, fermi_hf, spa_spectrum.second)) << endl;
    // }
    
  // }
  
}