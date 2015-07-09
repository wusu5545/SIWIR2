#include "functions.h"

extern boost::shared_ptr<double[]> minxyz;
extern boost::shared_ptr<double[]> maxxyz;
extern double r_cut;
extern double epsilon;
extern double sigma;

//data file
extern vector<double> x;
extern boost::shared_ptr<double[]> F;
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
	  boost::shared_ptr<int[]> cell_neighbor = boost::make_shared<int[]>(3,0);
	  //offset the distance between the particles
	  boost::shared_ptr<double[]> offset = boost::make_shared<double[]>(3,0);
	  cell_neighbor[0] = cell_num(i,0) + k;
	  cell_neighbor[1] = cell_num(i,1) + l;
	  cell_neighbor[2] = cell_num(i,2) + m;
	  for (size_t j;j<3;++j){
	    if (cell_neighbor[j]<0){
	      cell_neighbor[j] += cells_of_xyz[j];
	      offset[j] = -boxsize[j];
	    }
	    if (cell_neighbor[j]>=int(cells_of_xyz[j])){
	      cell_neighbor[j] -= cells_of_xyz[j];
	      offset[j] = boxsize[j];
	    }
	  }
	  int cell_neighbor_num = cell_neighbor[0]+cell_neighbor[1]*cells_of_xyz[0]+
				  cell_neighbor[2]*cells_of_xyz[0]*cells_of_xyz[1];
	  int p = cells[cell_neighbor_num];
	  //p = -1, no more particles in the cell
	  while (p!=-1){
	    if (p!=(int)i){//if i is not the particles in qusetion
	      boost::shared_ptr<double[]> r_ij = boost::make_shared<double[]>(3,0.0);
	      double r_ij_norm = 0;
	      for (size_t j=0;j<3;++j){
		r_ij[j] = x(i,j) - (x(p,j) + offset[j]);
		r_ij_norm += r_ij[j] * r_ij[j];
	      }
	      if (sqrt(r_ij_norm) <= r_cut)
		for (size_t j=0;j<3;++j){
		  double pow6 = (sigma*sigma/r_ij_norm)*(sigma*sigma/r_ij_norm)*(sigma*sigma/r_ij_norm);
		  F(i,j) += 24.0*epsilon/r_ij_norm*pow6*(2*pow6-1)*r_ij[j];
		}
	    }
	    p = particles[p];
	  }
	}
  }    
}
