#include "functions.h"

extern double sigma;

vector<double> InversePowerItration(const vector<Vertex> & A,double &eigenvalue)
{
  vector<double> eigenvector(A.size(),0.0);
  matrixvector M;
  eigenvector[0] = 1.0;
  eigenvalue = 0;
  M = A*eigenvector;
  for (size_t i = 0;i<100000000;++i)
  {
    eigenvector = ConjGrad(A,M.mass);
    eigenvector = 1.0/norm(eigenvector)*eigenvector;
    M = A*eigenvector;
    eigenvalue = eigenvector*M.stiff/(eigenvector*M.mass);
    if (norm(M.stiff-eigenvalue*M.mass)<sigma)
      break;
  }
  return eigenvector;
}