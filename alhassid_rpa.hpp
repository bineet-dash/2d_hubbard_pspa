#ifndef _ALHASSID_RPA_HPP_INCLUDED
#define _ALHASSID_RPA_HPP_INCLUDED

#include "pspa.hpp"
#include <chrono>


using namespace std::chrono;
typedef Matrix <cd, Dynamic, Dynamic, RowMajor> RMatrixXcd;

/* double get_pspa_F(MatrixXd u, VectorXd hf, double T)
{
	vector <MatrixXd> vt;
	for(int it=0; it<L; it++)
	{
		MatrixXd v_i = MatrixXd::Zero(2*L,2*L);
		v_i(it,it) = 1; v_i(it+L, it+L) = -1;
		MatrixXd v_i_transformed = u.adjoint()*v_i*u;
		vt.push_back(v_i_transformed);
	}

	double mu = get_mu(T, hf);
	VectorXd fermi_hf = VectorXd::Zero(hf.size());
	for(int it=0; it<hf.size(); it++)
	{
		fermi_hf(it) = fermi_fn(hf(it)-mu,T);
	}

	int r_max = int(abs( (hf(hf.size()-1)-hf(0))/T )) ; //omega_max = (2r_max+1)*pi*T= \delta_ij_max
	double final_det_r = 0;
	// ofstream fout("r_seqMst.dat");

	for(int matsubara_r = 0; matsubara_r < 3*r_max; matsubara_r++)
	{
		double omega_r = (2* matsubara_r +1)*M_PI*T;
		MatrixXcd rpa = MatrixXcd::Identity(L,L);
		for(int alpha=0; alpha<L; alpha++)
		{
			for(int alpha_prime=0; alpha_prime<L; alpha_prime++)
			{
				for(int i=0; i<hf.size(); i++)
				{
					for(int j=0; j<hf.size(); j++)
					{
						cd num = U/2*(vt.at(alpha_prime))(j,i)*(vt.at(alpha))(i,j)*(fermi_hf(i)-fermi_hf(j));
						cd denom = cd(hf(i)-hf(j),omega_r);
						rpa(alpha,alpha_prime) += num/denom;
					}
				}
			}
		}
		final_det_r += log( real(rpa.determinant()) );
		// fout << matsubara_r << " " << -log(real(rpa.determinant())) << endl;
	}

	return T*final_det_r;
} */

#endif