#ifndef _PSPA_HPP_INCLUDED
#define _PSPA_HPP_INCLUDED

#include <iostream>
#include <cmath>
#include <fstream>
#include "lib/lapacke.h"
#include <vector>
#include <Eigen/Dense>
#include <iomanip>
#include <string>
#include <chrono>

using namespace std;
using namespace Eigen;
using namespace std::chrono;

typedef std::complex <double> cd;
typedef pair <double,double> pdd;

extern double t;
extern double U;
extern int size;

double t=1;
double U;
int size;

inline double del(int a1, int a2){return (a1==a2)?1:0;}
inline cd jn(cd z){return conj(z);}
inline double Sqr(double x){return x*x;}
inline cd filter_cd(cd x){return (abs(x)<1e-4)?0.0:x;}
inline double filter_d(double x) {return (abs(x)<1e-4)?0.0:x;}
inline double filter_tol_d(double x, double tolerance=1e-4) {return (abs(x)<tolerance)?0.0:x;}
inline double fermi_fn(double e_minus_mu, double T) {return (std::isinf(exp(e_minus_mu/T)))? 0: 1/(exp(e_minus_mu/T)+1);}
inline int xc(int i, int sigma=1){return (sigma==1)?floor(i/size):floor((i-size*size)/size);}
inline int yc(int j, int sigma=1){return (sigma==1)?j%size:(j-size*size)%size;}
inline int index(int x, int y, int sigma=1){return (sigma==1)?x*size+y:x*size+y+size*size;}
inline int periodic(int a, int b, int lim) {int rem = (a+b)%lim; if(rem>=0) return rem; else return rem+lim;}

#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define MASK 123459876
double ran0(long *idum)
{
   long  k;
   double ans;

   *idum ^= MASK;
   k = (*idum)/IQ;
   *idum = IA*(*idum - k*IQ) - IR*k;
   if(*idum < 0) *idum += IM;
   ans=AM*(*idum);
   *idum ^= MASK;
   return ans;
}
#undef IA
#undef IM
#undef AM
#undef IQ
#undef IR
#undef MASK

bool zheev_cpp(MatrixXcd& A, vector<double>& lambda, char eigenvec_choice='N')
{
  int N = A.cols();
  int LDA = A.outerStride();
  int INFO = 0;
  double* w = new  double [N];
  char Nchar = 'N';
  char Vchar = 'V';
  char Uchar = 'U';
  int LWORK = int(A.size())*4;
  __complex__ double* WORK= new __complex__ double [LWORK];
  double* RWORK = new double [3*LDA];

  zheev_( &eigenvec_choice, &Uchar, &N, reinterpret_cast <__complex__ double*> (A.data()), &LDA, w, WORK, &LWORK, RWORK, &INFO );

  lambda.clear();
  for(int i=0; i<N; i++) lambda.push_back(w[i]);

  delete[] w; delete[] RWORK; delete[] WORK;
  return INFO==0;
}

bool zgeev_cpp(MatrixXcd& A, vector<double>& lambda, char eigenvec_choice='N')
{  
  int N = A.cols();
  int LDA = A.outerStride();
  int INFO = 0;
  __complex__ double* w = new __complex__ double [N];
  __complex__ double* vl;
  __complex__ double* vr;
  char Nchar = 'N';
  char Vchar = 'V';
  char Uchar = 'U';
  int LWORK = pow(2, N);
  __complex__ double* WORK= new __complex__ double [LWORK];
  double* RWORK = new double [LWORK];
  
  zgeev_(&Nchar, &eigenvec_choice, &N, reinterpret_cast <__complex__ double*> (A.data()), &LDA, w, vl, &LDA, vr, &LDA, WORK, &LWORK, RWORK, &INFO );

  lambda.clear();
  for(int i=0; i<N; i++) lambda.push_back(__real__ w[i]);

  delete[] w; delete[] RWORK; delete[] WORK;
  return INFO==0;
}

vector <double> stdEigenvalues(MatrixXcd A, bool (*diagonalization_routine)(MatrixXcd&, vector <double>&, char)=&zheev_cpp)
{
  std::vector<double> lambda; 
  if(diagonalization_routine(A,lambda,'N')) return lambda;
}

VectorXd Eigenvalues(MatrixXcd A, bool (*diagonalization_routine)(MatrixXcd&, vector <double>&, char)=&zheev_cpp)
{
  std::vector<double> lambda;
  if(diagonalization_routine(A,lambda,'N'))
 	{
		Map<ArrayXd> b(lambda.data(),lambda.size());
  	return b;
	}
}

MatrixXcd Eigenvectors(MatrixXcd A, bool (*diagonalization_routine)(MatrixXcd&, vector <double>&, char)=&zheev_cpp)
{
	std::vector<double> lambda;
  if(diagonalization_routine(A,lambda,'V')) return A; 
}

pair<MatrixXcd, vector<double>> stdEigenspectrum(MatrixXcd A, bool (*diagonalization_routine)(MatrixXcd&, vector <double>&, char)=&zheev_cpp)
{
  std::vector<double> lambda;
  if(diagonalization_routine(A,lambda,'V')) return make_pair(A,lambda);
}

pair<MatrixXcd, VectorXd> Eigenspectrum(MatrixXcd A, bool (*diagonalization_routine)(MatrixXcd&, vector <double>&, char)=&zheev_cpp)
{
  std::vector<double> lambda;
  if(diagonalization_routine(A,lambda,'V'))
 	{
    Map<ArrayXd> b(lambda.data(),lambda.size());
    return make_pair(A,b);
	}
}

void ising_sigma_generate(MatrixXd& suggested_randsigma, int lattice_index, long & idum)
{
  if(ran0(&idum)<=0.5) suggested_randsigma(lattice_index,2) *= -1;
}

MatrixXcd construct_h0_2d(void)
{
  MatrixXcd Mc = MatrixXcd::Zero(2*size*size, 2*size*size);
  for(int i=0; i<size*size; i++)
  {
    for(int j=0; j<size*size; j++)
    {
      int ix = xc(i); int iy = yc(i); int jx = xc(j); int jy = yc(j);
      if((ix == periodic(jx,1,size)|| ix == periodic(jx,-1,size)) && iy == jy ) Mc(i,j) = cd(-t,0);
			if(ix == jx && (iy == periodic(jy,1,size) || iy == periodic(jy,-1,size))) Mc(i,j) = cd(-t,0);				
    }
  }
  for(int i=size*size; i<2*size*size; i++)
  {
    for(int j=size*size; j<2*size*size; j++)
    {
      int ix = xc(i,-1); int iy = yc(i,-1); int jx = xc(j,-1); int jy = yc(j,-1);
      if((ix == periodic(jx,1,size)|| ix == periodic(jx,-1,size)) && iy == jy ) Mc(i,j) = cd(-t,0);
			if(ix == jx && (iy == periodic(jy,1,size) || iy == periodic(jy,-1,size))) Mc(i,j) = cd(-t,0);				
    }
  }
  return Mc;
}

MatrixXcd matrixelement_sigmax_2d(MatrixXd randsigma)
{
  MatrixXcd Mcx = MatrixXcd::Zero(2*size*size,2*size*size);
  for(int row=0; row<size*size; row++)
  {
    Mcx(row,row+size*size) = cd(randsigma(row,0),0);
    Mcx(row+size*size,row) = cd(randsigma(row,0),0);
  }
  return Mcx;
}

MatrixXcd matrixelement_sigmay_2d(MatrixXd randsigma)
{
  MatrixXcd Mcy = MatrixXcd::Zero(2*size*size,2*size*size);
  for(int row=0; row<size*size; row++)
  {
    Mcy(row,row+size*size) = cd(0,-randsigma(row,1));
    Mcy(row+size*size,row) = cd(0,randsigma(row,1));
  }
  return Mcy;
}

MatrixXcd matrixelement_sigmaz_2d(MatrixXd randsigma)
{
  MatrixXcd Mcz = MatrixXcd::Zero(2*size*size,2*size*size);
  for(int row=0; row<size*size; row++)  Mcz(row,row)= cd(randsigma(row,2),0);
  for(int row=size*size; row<2*size*size; row++)  Mcz(row, row)=cd(-randsigma(row-size*size,2),0);
  return Mcz;
}


double get_spi(MatrixXd sigma )
{
  double sq = 0.0;
  for(int i=0; i<size*size; i++)
  {
    for(int j=0; j<size*size; j++)
    {
      sq += sigma(i,2)*sigma(j,2)*pow(-1,xc(i)-xc(j))*pow(-1,yc(i)-yc(j)) /pow(size,4);
    } 
  }
  return sq;
}

double get_mu(double temperature, std::vector<double> v)
{
  sort (v.begin(), v.end());
  double bisection_up_lim = v.back();
  double bisection_low_lim = v.front();

  double mu, no_of_electrons; int count=0;
  double epsilon = 0.000001;

  for(; ;)
  {
    no_of_electrons=0;
    mu = 0.5*(bisection_low_lim+bisection_up_lim);

    for(auto it = v.begin(); it!= v.end(); it++)
    {
      double fermi_func = 1/(exp((*it-mu)/temperature)+1);
      no_of_electrons += fermi_func;
    }
    if(abs(no_of_electrons-size*size) < epsilon)
    {
      return mu; break;
    }
    else if(no_of_electrons > size*size+epsilon)
    {
       if(abs(bisection_up_lim-v.front())<0.001){return mu; break;}
       else {bisection_up_lim=mu;}
    }
    else if(no_of_electrons < size*size-epsilon)
    {bisection_low_lim=mu;}
  }
}

double get_mu(double temperature, VectorXd v)
{
  vector<double> stdv (v.data(),v.data()+v.size());
  return get_mu(temperature, stdv);
}

double spa_free_energy(Eigen::VectorXd spa_eivals, double T)
{
  double mu = get_mu(T, spa_eivals);
  double beta = 1/T;
  double spa_F = 0.0;
  for(int i=0; i<spa_eivals.size(); i++)
  {
    spa_F += (-beta*(spa_eivals(i)-mu) > 4.0)?  (spa_eivals(i)-mu):-T*log(1+exp(-beta*(spa_eivals(i)-mu)));
  }
  return spa_F+mu*size*size;
}

/* double spa_internal_energy(MatrixXcd Mc, double temperature)
{
  std::vector<double> eigenvalues;
  zheev_cpp(Mc, eigenvalues, 'N');
  sort(eigenvalues.begin(),eigenvalues.end());
  double mu = get_mu(temperature, eigenvalues);

  double internal_energy=0.0 ; double e_min = eigenvalues.front();
  for(auto it=eigenvalues.begin(); it!= eigenvalues.end(); it++)
  {
    internal_energy += (*it)/(exp((*it-mu)/temperature)+1);
  }
  return internal_energy/(size*size);
} */

double get_pspa_F(MatrixXd u, VectorXd hf, double T)
{
	vector <MatrixXd> vt;
	for(int it=0; it<size*size; it++)
	{
		MatrixXd v_i = MatrixXd::Zero(2*size*size,2*size*size);
		v_i(it,it) = 1; v_i(it+size*size, it+size*size) = -1;
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
		MatrixXcd rpa = MatrixXcd::Identity(size*size,size*size);

		milliseconds begin_ms, end_ms;
		begin_ms = duration_cast< milliseconds >(system_clock::now().time_since_epoch());
		for(int alpha=0; alpha<size*size; alpha++)
		{
			for(int alpha_prime=0; alpha_prime<size*size; alpha_prime++)
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
		end_ms = duration_cast< milliseconds >(system_clock::now().time_since_epoch());
		final_det_r += log( real(rpa.determinant()) );
	}

	return T*final_det_r;
}


pair <double, double> get_spa_pspa_F(MatrixXd u, VectorXd spa_eivals, double temperature)
{
	double spa_F =  spa_free_energy(spa_eivals, temperature);
	double pspa_F = spa_F + get_pspa_F(u, spa_eivals, temperature)/(size*size);
	return make_pair(spa_F, pspa_F);
}

VectorXd inttobin(int theValue)
{
  VectorXd v(size*size);
  for (int i = 0; i < size*size; ++i)  v(size*size-1-i) = theValue & (1 << i) ? 1 : 0;
  return v;
}

VectorXd get_field(int i)
{
  VectorXd raw = inttobin(i);
  for(int i=0; i<raw.size(); i++) raw(i) = (raw(i)==0)?-1:1;
  return raw;
}

#endif
