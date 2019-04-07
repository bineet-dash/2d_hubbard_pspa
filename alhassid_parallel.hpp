#ifndef _ALHASSID_PARALLEL_HPP_INCLUDED_
#define _ALHASSID_PARALLEL_HPP_INCLUDED_

#include "pspa.hpp"
#include <chrono>
#include <mpi.h>


using namespace std::chrono;
typedef Matrix <cd, Dynamic, Dynamic, RowMajor> RMatrixXcd;

int NO_MC_TRIALS = 300;


double f_det(int matsubara_r, double temperature, const VectorXd& spa_eivals, const VectorXd& fermi_hf, const vector<MatrixXd>& vt)
{
	int L = size*size;
	double omega_r = (2* matsubara_r +1)*M_PI*temperature;
	MatrixXcd rpa = MatrixXcd::Identity(L,L);

	for(int alpha=0; alpha<L; alpha++)
	{
		for(int alpha_prime=0; alpha_prime<L; alpha_prime++)
		{
			for(int i=0; i<spa_eivals.size(); i++)
			{
				for(int j=0; j<spa_eivals.size(); j++)
				{
					cd num = U/2*(vt.at(alpha_prime))(j,i)*(vt.at(alpha))(i,j)*(fermi_hf(i)-fermi_hf(j));
					cd denom = cd(spa_eivals(i)-spa_eivals(j),omega_r);
					rpa(alpha,alpha_prime) += num/denom;
				}
			}
		}
	}
	return 1/rpa.real().determinant();
}

double profile_f_det(int matsubara_r, double temperature, const VectorXd& spa_eivals, const VectorXd& fermi_hf, const vector<MatrixXd>& vt, int pRank)
{
	int L = size*size;
	double omega_r = (2* matsubara_r +1)*M_PI*temperature;
	MatrixXcd rpa = MatrixXcd::Identity(L,L);
	milliseconds begin_ms, end_ms;
   begin_ms = duration_cast< milliseconds >(system_clock::now().time_since_epoch());

	for(int alpha=0; alpha<L; alpha++)
	{
	  begin_ms = duration_cast< milliseconds >(system_clock::now().time_since_epoch());
		for(int alpha_prime=0; alpha_prime<L; alpha_prime++)
		{
			for(int i=0; i<spa_eivals.size(); i++)
			{
				for(int j=0; j<spa_eivals.size(); j++)
				{
					cd num = U/2*(vt.at(alpha_prime))(j,i)*(vt.at(alpha))(i,j)*(fermi_hf(i)-fermi_hf(j));
					cd denom = cd(spa_eivals(i)-spa_eivals(j),omega_r);
					rpa(alpha,alpha_prime) += num/denom;
				}
			}
		}
  	end_ms = duration_cast< milliseconds >(system_clock::now().time_since_epoch());	
		if(pRank==0) cout << "row " << alpha << " takes= " << double((end_ms-begin_ms).count())/1000 << endl;
	}
	end_ms = duration_cast< milliseconds >(system_clock::now().time_since_epoch());	

	return  1/rpa.real().determinant();
}

double det_sum_MC(double temperature, int r_max, const VectorXd& spa_eivals, const VectorXd& fermi_hf,  const vector <MatrixXd>& vt)
{
	int num_procs, my_rank;
	MPI_Comm_size (MPI_COMM_WORLD, &num_procs);
	MPI_Comm_rank (MPI_COMM_WORLD, &my_rank);
	
	srand(time(NULL)+my_rank);
	long idum = time(NULL)+my_rank;
	int selected_MC = 0;

	double fr_max = log(f_det(0,temperature, spa_eivals, fermi_hf, vt));

	int local_selected_MC = 0;
	for(int i=0; i<NO_MC_TRIALS/num_procs; i++)
	{
		int r_rand = rand()%r_max;
		double exp_fr_rand = f_det(r_rand, temperature, spa_eivals, fermi_hf, vt);
		if(exp_fr_rand < 0) continue;
		double y_rand = ran0(&idum)*fr_max;
		if(y_rand <= log(exp_fr_rand))
		{
			local_selected_MC ++;
		}
	}	
	MPI_Reduce(&local_selected_MC, &selected_MC,1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

	if(my_rank==0)
	{
		double MC_total_det_r = r_max*fr_max*(double(selected_MC)/double(NO_MC_TRIALS));
		return -temperature*MC_total_det_r;
	}
}

double det_sum_parallel(double temperature, int r_max, const VectorXd& spa_eivals, const VectorXd& fermi_hf,  const vector <MatrixXd>& vt)
{
	double final_det_r = 0;

	int num_procs, my_rank;
	MPI_Comm_size (MPI_COMM_WORLD, &num_procs);
	MPI_Comm_rank (MPI_COMM_WORLD, &my_rank);

	int up_limit = 3*r_max + (num_procs - 3*r_max%num_procs);

	if(up_limit% num_procs!=0) 
	{
		cout << up_limit << " " << num_procs << up_limit%num_procs << endl;
		exit(32);
	}

	int elem_per_proc = up_limit/num_procs;
	double local_det_r = 0.0;

	for(int it = 0; it < elem_per_proc; it++)
	{
		int matsubara_r = my_rank*elem_per_proc+it;
		double f_r = f_det(matsubara_r, temperature, spa_eivals, fermi_hf, vt);
		if(f_r > 1e-4)
		{
			local_det_r += log( f_r );
		}
		else
		{
			cerr << "Negative f_r detected. " << endl;
		}
		
	}
	MPI_Comm_set_errhandler(MPI_COMM_WORLD, MPI_ERRORS_RETURN);
	int ierr = MPI_Reduce(&local_det_r, &final_det_r,1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

	if(my_rank==0)
	{
		return -temperature*final_det_r;
	}

}

double det_sum_explicit(double temperature, int r_max, const VectorXd& spa_eivals, const VectorXd& fermi_hf,  const vector <MatrixXd>& vt)
{
	double final_det_r = 0;
	for(int matsubara_r = 0; matsubara_r < 3*r_max; matsubara_r++)
	{
		double f_r = f_det(matsubara_r, temperature, spa_eivals, fermi_hf, vt);
		if(f_r>0)
		{
		 final_det_r += log( f_r );
		}
		else
		{
			cerr << "Negative f_r detected. " << endl;
		}
	}
	return -temperature*final_det_r;
}

double pspa_free_energy(double temperature, const VectorXd& spa_eivals, const MatrixXd& u, double fill=1.0)
{
	int r_max = int(abs( (spa_eivals(spa_eivals.size()-1)-spa_eivals(0))/temperature )) ; //omega_max = (2r_max+1)*pi*T= \delta_ij_max
	int L = size*size;
	
	vector <MatrixXd> vt;
	for(int it=0; it<L; it++)
	{
		MatrixXd v_i = MatrixXd::Zero(2*L,2*L);
		v_i(it,it) = 1; v_i(it+L, it+L) = -1;
		MatrixXd v_i_transformed = u.adjoint()*v_i*u;
		vt.push_back(v_i_transformed);
	}

	double mu = get_mu(temperature, spa_eivals, fill);
	VectorXd fermi_hf = VectorXd::Zero(spa_eivals.size());
	for(int it=0; it< spa_eivals.size(); it++)
	{
		fermi_hf(it) = fermi_fn(spa_eivals(it)-mu, temperature);
	}

	int pRank, num_procs;
	MPI_Comm_size (MPI_COMM_WORLD, &num_procs);
	MPI_Comm_rank (MPI_COMM_WORLD, &pRank);

	double pspa_F = 0; 
	if(r_max < num_procs)
	{
		pspa_F = det_sum_explicit(temperature, r_max, spa_eivals, fermi_hf, vt);
	}
	else if(r_max > num_procs && r_max < int(NO_MC_TRIALS/3))
	{
		pspa_F = det_sum_parallel(temperature, r_max, spa_eivals, fermi_hf, vt);
	}
	else
	{
		pspa_F = det_sum_MC(temperature, r_max, spa_eivals, fermi_hf, vt);
	}
	return pspa_F;
}



double profile_pspa_free_energy(double temperature, const VectorXd& spa_eivals, const MatrixXd& u, int pRank, double fill=1.0)
{
   milliseconds begin_ms, end_ms;
	int r_max = int(abs( (spa_eivals(spa_eivals.size()-1)-spa_eivals(0))/temperature )) ; //omega_max = (2r_max+1)*pi*T= \delta_ij_max
	int L = size*size;
	
	if(pRank==0) cout << "rmax = " << r_max << endl;
	
	vector <MatrixXd> vt;
	for(int it=0; it<L; it++)
	{
		MatrixXd v_i = MatrixXd::Zero(2*L,2*L);
		v_i(it,it) = 1; v_i(it+L, it+L) = -1;
		MatrixXd v_i_transformed = u.adjoint()*v_i*u;
		vt.push_back(v_i_transformed);
	}

	double mu = get_mu(temperature, spa_eivals, fill);
	VectorXd fermi_hf = VectorXd::Zero(spa_eivals.size());
	for(int it=0; it< spa_eivals.size(); it++)
	{
		fermi_hf(it) = fermi_fn(spa_eivals(it)-mu, temperature);
	}

	double pspa_F = 0; 
   begin_ms = duration_cast< milliseconds >(system_clock::now().time_since_epoch());
	// profile_f_det(3, temperature, spa_eivals, fermi_hf, vt, pRank);

	// if(r_max < 32)
	// {
	// 	pspa_F = det_sum_explicit(temperature, r_max, spa_eivals, fermi_hf, vt);
	// 	if(pRank==0){cout << "Took " <<  double((end_ms-begin_ms).count())/1000.0 << " seconds. " << "\r"; cout.flush();}
	// }
	if(/* r_max > 32 && */ r_max < 100)
	{
		pspa_F = det_sum_parallel(temperature, r_max, spa_eivals, fermi_hf, vt);
	}
	else
	{
		pspa_F = det_sum_MC(temperature, r_max, spa_eivals, fermi_hf, vt);
	}

	end_ms = duration_cast< milliseconds >(system_clock::now().time_since_epoch());
	if(pRank==0) cout << "Took " <<  double((end_ms-begin_ms).count())/1000.0 << " seconds. " << endl; 
	return pspa_F;
}

#endif