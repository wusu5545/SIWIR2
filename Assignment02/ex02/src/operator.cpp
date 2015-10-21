#include "functions.h"

matrixvector operator*(const vector<Vertex> &A,const vector<double> &u)
{
  matrixvector res;
  double stiffness,mass;
  for (size_t i=0;i<A.size();++i){
    stiffness=0.0;
    mass=0.0;
    for (auto itr=A[i].b();itr!=A[i].e();++itr){
      stiffness += u[itr->first]*itr->second.stiffness;
      mass += u[itr->first]*itr->second.mass;
    }
    res.stiff.push_back(stiffness);
    res.mass.push_back(mass);
  }
  return res;
}

vector<double> operator*(const double & a,const vector<double> &b)
{
  vector<double> res;
  for (size_t i=0;i<b.size();++i)
    res.push_back(a*b[i]);
  return res;
}

double operator*(const vector<double> &a,const vector<double> &b)
{
  double res=0.0;
  for(size_t i=0;i<b.size();++i)
    res += a[i]*b[i];
  return res;
}

vector<double> operator+(const vector<double> &a,const vector<double> &b)
{
  vector<double> res;
  for(size_t i=0;i<b.size();++i)
    res.push_back(a[i]+b[i]);
  return res;
}

vector<double> operator-(const vector<double> &a,const vector<double> &b)
{
  vector<double> res;
  for(size_t i=0;i<b.size();++i)
    res.push_back(a[i]-b[i]);
  return res;
}