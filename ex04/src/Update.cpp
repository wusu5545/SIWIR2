#include "functions.h"

extern double delta_t;
extern vector<double> x;
extern vector<double> v;
extern boost::shared_ptr<double[]> F;
extern boost::shared_ptr<double[]> F_old;
extern vector<double> mass;
extern size_t N;

void Position_update()
{
  for (size_t i = 0;i<N;++i)
    for (size_t j = 0;j<3;++j)
      x(i,j) += delta_t*v(i,j) + F(i,j)*delta_t*delta_t/(2*mass[i]);
}

void Velocity_update()
{
  for (size_t i = 0;i<N;++i)
    for (size_t j = 0;j<3;++j)
      v(i,j) += (F_old(i,j) + F(i,j))*delta_t/(2*mass[i]);
}