#include<iostream>
#include<fstream>
#include<stdlib.h>
#include<boost/make_shared.hpp>
#include<boost/shared_ptr.hpp>
#include<string>
#include<iomanip> 
#include<map>
#include<vector>
#include"vertex.h"

using namespace std;
using namespace boost;
using std::vector;

matrixvector operator*(const vector<Vertex> &A,const vector<double> &x);
vector<double> operator*(const double & a,const vector<double> &b);
double operator*(const vector<double> &a,const vector<double> &b);
vector<double> operator+(const vector<double> &a,const vector<double> &b);
vector<double> operator-(const vector<double> &a,const vector<double> &b);

double norm(const vector<double> &v);
vector<double> ConjGrad(const vector<Vertex> &A,const vector<double> & b);
vector<double> InversePowerItration(const vector<Vertex> & A,double &eigenvalue);
void Waveguide( const string filename);
double ksq(double x,double y);