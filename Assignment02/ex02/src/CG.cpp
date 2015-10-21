#include "functions.h"

vector<double> ConjGrad(const vector<Vertex> &A,const vector<double> & b)
{
  vector<double> r;
  vector<double> x(A.size(),0.0);
  vector<double> p;
  vector<double> Ap;
  double alpha;
  double rsnew;
  
  r = b - (A*x).stiff;
  p = r;
  double rsold = r*r;
  
  for (size_t i=0;i<100000000;++i)
  {
    Ap = (A*p).stiff;
    alpha = rsold / (p*Ap);
    x = x + alpha*p;
    r = r - alpha*Ap;
    rsnew = r*r;
    if (sqrt(rsnew)<1e-14)
      break;
    p = r + rsnew/rsold * p;
    rsold = rsnew;
  }
  
  return x;
}