#include "functions.h"

double norm(const vector<double> &v)
{
  double res = 0.0;
  for (size_t i= 0;i<v.size()-1;++i)
    res += v[i]*v[i];
  res = sqrt(res);
  return res;
}