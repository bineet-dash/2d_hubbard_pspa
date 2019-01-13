#ifndef _SPA_LIBRARY_HPP_INCLUDED_
#define _SPA_LIBRARY_HPP_INCLUDED_

#include <iostream>
#include <cmath>
#include <fstream>
#include <Eigen/Dense>
#include <complex>
#include <iomanip>
#include <string>
#include <chrono>
#include <lapacke.h>

using namespace std;
using namespace Eigen;
using namespace std::literals;
using namespace std::chrono;


typedef std::complex<double> cd;

extern int size;
extern double t;
extern double U;


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

double t=1;
double U;
int size;

double Sqr(cd z){return norm(z);}
double filter(double x) {if(abs(x)<1e-7) return 0.0; else return x;}
void filter(std::vector<double>& v) {for(int i=0; i<v.size(); i++)  v[i]=filter(v[i]); }
void filter(VectorXd& v) {for(int i=0; i<v.size(); i++)  v(i)=filter(v(i));}
inline int xc(int i, int sigma=1){return (sigma==1)?floor(i/size):floor((i-size*size)/size);}
inline int yc(int j, int sigma=1){return (sigma==1)?j%size:(j-size*size)%size;}
inline int index(int x, int y, int sigma=1){return (sigma==1)?x*size+y:x*size+y+size*size;}
inline int periodic(int a, int b, int lim) {int rem = (a+b)%lim; if(rem>=0) return rem; else return rem+lim;}


bool diagonalize(MatrixXcd A, vector<double>& lambda)
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
  zheev_( &Nchar, &Uchar, &N, reinterpret_cast <__complex__ double*> (A.data()), &LDA, w, WORK, &LWORK, RWORK, &INFO );

  lambda.clear();
  for(int i=0; i<N; i++) lambda.push_back(w[i]);

  delete[] w; delete[] RWORK; delete[] WORK;
  return INFO==0;
}

VectorXd Eigenvalues(MatrixXcd A)
{
  std::vector<double> lambda;
  bool ret_value = diagonalize(A,lambda);
  Map<ArrayXd> b(lambda.data(),lambda.size());
  return b;
}

void ising_sigma_generate(MatrixXd& suggested_randsigma, int lattice_index, long & idum)
{
  if(ran0(&idum)<=0.5) suggested_randsigma(lattice_index,2) *= -1;
}

/* void sigma_generate(MatrixXd& randsigma, int i, long & idum, double T)
{
  double radius, u, theta;
  // radius = (T<0.5)? 0.5+ran0(&idum): ( min(2*T,8.0)*ran0(&idum) );
  radius = 0.5+ran0(&idum);
  u = 2*ran0(&idum)-1;
  theta  = 2*3.1416*ran0(&idum);
  randsigma(i,0)= radius*sqrt(1-pow(u,2))*cos(theta); 
  randsigma(i,1)= radius*sqrt(1-pow(u,2))*sin(theta);
  randsigma(i,2)= radius*u;
}
*/

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

double spa_free_energy(MatrixXcd Mc, double temperature)
{
  std::vector<double> eigenvalues;
  diagonalize(Mc, eigenvalues);
  sort(eigenvalues.begin(),eigenvalues.end());
    // for(auto it=eigenvalues.begin(); it!=eigenvalues.end(); it++) cout << filter(*it) << " "; cout << endl;

  double free_energy = 0; double ekt =0;
  double mu = get_mu(temperature, eigenvalues);
    // cout << "mu= " << mu << endl;

  for(auto it=eigenvalues.begin(); it!= eigenvalues.end(); it++)
  {
    ekt = (*it-mu)/temperature;
    if(!isinf(exp(-ekt))) free_energy += -temperature*log(1+exp(-ekt));
    else  free_energy += (*it-mu);
  }
  //  cout << free_energy << endl;
  return free_energy+size*size*mu;
}

double spa_internal_energy(MatrixXcd Mc, double temperature)
{
  std::vector<double> eigenvalues;
  diagonalize(Mc, eigenvalues);
  sort(eigenvalues.begin(),eigenvalues.end());

  double internal_energy = 0; double ekt =0;
  double mu = get_mu(temperature, eigenvalues);
    // cout << "mu= " << mu << endl;

  for(auto it=eigenvalues.begin(); it!= eigenvalues.end(); it++)
  {
    internal_energy += *it/(exp((*it-mu)/temperature)+1);
  }
  //  cout << free_energy << endl;
  return internal_energy;
}

#endif
