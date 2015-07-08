#include "functions.h"

extern boost::shared_ptr<double[]> minxyz;
extern boost::shared_ptr<double[]> maxxyz;
extern double r_cut;
extern double epsilon;
extern double sigma;

//data file
extern vector<double> x;
extern boost::shared_ptr<double[]> F;
extern boost::shared_ptr<double[]> F_old;
extern boost::shared_ptr<double[]> r;
extern vector<double> mass;

//derived parameters
extern size_t N;
extern boost::shared_ptr<double[]> boxsize;
extern boost::shared_ptr<size_t[]> cells_of_xyz;
extern boost::shared_ptr<double[]> cellsize;


//link cell
extern boost::shared_ptr<int[]> cells;
extern boost::shared_ptr<int[]> particles;
extern boost::shared_ptr<size_t[]> cell_num;

void ForceCalculation()
{
  for (size_t i=0;i<N;++i){
    for (size_t j=0;j<3;++j)
      F(i,j) = 0;
    for (int k=-1;k<2;++k)
      for (int l=-1;l<2;++l)
	for (int m=-1;m<2;++m){
	  
	}
  }
      
}
