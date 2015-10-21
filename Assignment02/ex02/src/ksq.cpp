#include<cmath>

extern double delta;

double ksq(double x,double y)
{
  return (100.0 + delta)*exp(-50.0*(x*x+y*y)) - 100.0;
}