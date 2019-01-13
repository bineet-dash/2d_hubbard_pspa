#ifndef _EXTRA_HPP_INCLUDED_
#define _EXTRA_HPP_INCLUDED_

#include <fstream>
#include <Eigen/Dense>
#include <sys/types.h>
#include <sys/stat.h>
#include <ctime>
#include <iostream>
#include <iomanip>
#include <cstdio>
#include <string>

using namespace std;
using namespace Eigen;

void show_eigenvalues(MatrixXcd H)
{
  std::vector<double> eigenvalues;
  zheev_cpp(H, eigenvalues);
  for(auto it=eigenvalues.begin(); it!= eigenvalues.end(); it++) cout << *it << " ";
  cout << endl << endl;
}

void matrix_output(MatrixXcd Mc, MatrixXcd Mcx, MatrixXcd Mcy, MatrixXcd Mcz)
{
  cout << setprecision(3) <<  Mc << endl << endl << Mcx << endl << endl << Mcy << endl << endl << Mcz << endl << endl;
}

void generate_scriptname(string& scriptname) {
    time_t rawtime;
    struct tm * timeinfo;
    char buffer[80];
    time (&rawtime);
    timeinfo = localtime(&rawtime);
    strftime(buffer,sizeof(buffer),"%S-%M-%I-%Y-%m-%d",timeinfo);
    string str(buffer);
    string base= "wolframscripts/script";
    string extension=".nb";
    scriptname= base+str+extension;
}

string current_time_str(void)
{
  time_t rawtime;
  struct tm * timeinfo;
  char buffer[80];
  time (&rawtime);
  timeinfo = localtime(&rawtime);
  strftime(buffer,sizeof(buffer),"%S-%M-%I-%Y-%m-%d",timeinfo);
  string str(buffer);
  return str;
}

void eigenvalues_Mathematica(MatrixXcf Mc, ofstream& fout, string scriptname)
{
  int size = Mc.rows()/2;
  fout.open(scriptname);
  //fout << "#!/usr/local/bin/WolframScript -linewise -script" << endl;
  fout << "Print[Eigenvalues[N[{";
  for(int i=0; i<2*size; i++)
  {
    fout << "{ ";
    for(int j=0; j<2*size; j++)
      {
        fout << Mc.real()(i,j) << "+ I (" << Mc.imag()(i,j);
         if(j==2*size-1) fout << ") ";
         else fout << "),";
      }
   if(i==2*size-1) fout << "} ";
   else fout << "},";
  }
  fout << "}]]];" << endl;
  fout.close();
}

void make_exec(string str)
{
   const char *scriptname = str.c_str();
   chmod(scriptname, S_IRWXU);
   cout << "To generate the eigenvalues using Mathematica for verification, run "<< str << endl;
}

void generate_wlscript(MatrixXcf Mc, ofstream& fout)
{
  int size = Mc.rows()/2;
  string scriptname;
  generate_scriptname(scriptname);
  eigenvalues_Mathematica(Mc, fout, scriptname);
  make_exec(scriptname);
}

void progress_percent_desc(double initial_temp, double final_temp, double temperature)
{
   int progress_percent;
   progress_percent=int((final_temp-temperature)*100/(final_temp-initial_temp));
   cout.flush();
   cout << "\r [ "<< progress_percent+1 << "% ] ";
   for(int i=0; i<progress_percent; i++)
   {
     cout << "#";
   }
}

void progress_percent_asc(double initial_temp, double final_temp, double temperature)
{
   int progress_percent;
   progress_percent=int((temperature-initial_temp)*100/(final_temp-initial_temp));
   cout.flush();
   cout << "\r [ "<< progress_percent+1 << "% ] ";
   for(int i=0; i<progress_percent; i++)
   {
     cout << "#";
   }
}

void spinarrangement_2d_Mathematica_output(MatrixXd M, ofstream& nb, double lattice_separation = 1.0)
{

  nb << "Show[Graphics3D[{" << endl;
  for(int i=0; i< size; i++)
  {
    for(int j=0; j< size; j++)
    {
      double lattice_x = lattice_separation*i;
			double lattice_y = lattice_separation*j; 
			double spin_x = M(index(i,j),0) + lattice_x;
			double spin_y = M(index(i,j),1) + lattice_y;
      double spin_z = M(index(i,j),2);

      nb << "Arrow[{{" << lattice_x << ", " << lattice_y << ",0}, {" << spin_x << ","  << spin_y << "," << spin_z <<"}}]";
			if(i!=size-1 || j!=size-1) nb << ",";
    }
  }
  nb <<"}] ]" << endl << endl;
}

void show_time(milliseconds begin_ms, milliseconds end_ms, string s)
{
   long int t = (end_ms.count()-begin_ms.count())/1000;
    if (t<=60)
    { cout <<  s << " took " << t << " seconds." << endl; }
    else if (t>60 && t<3600)
    {
      int minutes=int(t/60); int seconds= t%minutes;
      cout << s << " took " << minutes << " minute and " << seconds << " seconds." << endl;
    }
    else if(t>=3600)
    {
      int hrs= int(t/3600); int minutes=int((t-3600*hrs)/60); int seconds= int(t-3600*hrs-60*minutes);
      cout << s << " took " << hrs << " hour, " << minutes << " minutes and " << seconds << " seconds. ";
    }
    else
    {cout << s << " took " << t << "time. Wrong t received.\n"; }
}

#endif
