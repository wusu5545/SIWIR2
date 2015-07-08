#include "functions.h"

//parameter file
extern boost::shared_ptr<double[]> minxyz;
extern boost::shared_ptr<double[]> maxxyz;

//data file
extern vector<double> x;

//derived parameters
extern size_t N;
extern boost::shared_ptr<double[]> boxsize;
extern boost::shared_ptr<size_t[]> cells_of_xyz;
extern boost::shared_ptr<double[]> cellsize;


//link cell
extern boost::shared_ptr<int[]> cells;
extern boost::shared_ptr<int[]> particles;
extern boost::shared_ptr<size_t[]> cell_num;

void LinkedCell()
{
  //treat boundary
  for (size_t i=0;i<N;++i)
    for (size_t j=0;j<3;++j){
      if (x(i,j)<minxyz[j])
	x(i,j) += boxsize[j];
      if (x(i,j)>=minxyz[j])
	x(i,j) -= boxsize[j];
    }
  
  for (size_t i=0;i<N;++i){
    for (size_t j=0;j<3;++j)
      cell_num(i,j) = x(i,j)/cellsize[j];
    size_t cell_num_total = cell_num(i,0)+cell_num(i,1)*cells_of_xyz[0]+
			    cell_num(i,2)*cells_of_xyz[1]*cells_of_xyz[0];
    swap(particles[i],cells[cell_num_total]);
  }
}
