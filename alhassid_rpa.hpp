#ifndef _ALHASSID_RPA_HPP_INCLUDED
#define _ALHASSID_RPA_HPP_INCLUDED

#include "rpa.hpp"
#include <chrono>


using namespace std::chrono;
typedef Matrix <cd, Dynamic, Dynamic, RowMajor> RMatrixXcd;

double rpa_det(double omega, vector <MatrixXd> vt, VectorXd hf, double T)
{
	MatrixXd rpa = MatrixXd::Identity(L,L);
	double mu = get_mu(T, hf);

	for(int alpha=0; alpha<L; alpha++)
	{
		for(int alpha_prime=0; alpha_prime<L; alpha_prime++)
		{
			for(int i=0; i<hf.size(); i++)
			{
				for(int j=0; j<hf.size(); j++)
				{
					rpa(alpha,alpha_prime) += U_prime/2*(vt.at(alpha_prime))(j,i)*(vt.at(alpha))(i,j)*(fermi_fn(hf(i)-mu,T)-fermi_fn(hf(j)-mu,T))/(hf(i)-hf(j)+omega);
				}
			}
		}
	}
	double result = rpa.determinant();
	return result;
}

double get_pspa_F(MatrixXd u, VectorXd hf, double T)
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
	
	for(int matsubara_r = 0; matsubara_r < 5*r_max; matsubara_r++)
	{
		double omega_r = (2* matsubara_r +1)*M_PI*T;
		MatrixXcd rpa = MatrixXcd::Identity(L,L);

		milliseconds begin_ms, end_ms;
		begin_ms = duration_cast< milliseconds >(system_clock::now().time_since_epoch());
		for(int alpha=0; alpha<L; alpha++)
		{
			for(int alpha_prime=0; alpha_prime<L; alpha_prime++)
			{
				for(int i=0; i<hf.size(); i++)
				{
					for(int j=0; j<hf.size(); j++)
					{
						cd num = U_prime/2*(vt.at(alpha_prime))(j,i)*(vt.at(alpha))(i,j)*(fermi_hf(i)-fermi_hf(j));
						cd denom = cd(hf(i)-hf(j),omega_r);
						rpa(alpha,alpha_prime) += num/denom;
					}
				}
			}
		}
		end_ms = duration_cast< milliseconds >(system_clock::now().time_since_epoch());
		final_det_r += log( real(rpa.determinant()) );
	}

	return T*final_det_r;
}



pair <double, double> get_spa_pspa_F(MatrixXd u, VectorXd spa_eivals, double temperature)
{
	double spa_F =  spa_free_energy(spa_eivals, temperature);
	double pspa_F = spa_F + get_pspa_F(u, spa_eivals, temperature)/L;
	return make_pair(spa_F, pspa_F);
}

VectorXd inttobin(int theValue)
{
  VectorXd v(L);
  for (int i = 0; i < L; ++i)  v(L-1-i) = theValue & (1 << i) ? 1 : 0;
  return v;
}

VectorXd get_field(int i)
{
  VectorXd raw = inttobin(i);
  for(int i=0; i<raw.size(); i++) raw(i) = (raw(i)==0)?-1:1;
  return raw;
}

#endif