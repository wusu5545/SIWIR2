#include <math.h>
#include  <iostream>
#include <fstream>
#include "FUNCTIONS.h"
#include <vector>
#include <iomanip>
#include <stdlib.h>

void solveMG(int l, std::vector<std::vector<double > >& u, std::vector<std::vector<double > >& f){
      
   //DEFINITION OF PARAMETERS
  int v1 = 3;
  int v2 = 2;
  int gamma =1;
  int miu = 1;
  std::size_t grid_size = pow(2,l)+1;

//   std::size_t grid_size_x = grid_size;
//   std::size_t grid_size_y = (grid_size - 1)/2 + 1;

//START OF SOLUTION

  for(int k=0; k<4; ++k){

      fmg(u, f, v1, v2, gamma, miu, l, grid_size);

	//PRINT OUT DESCRETE L2 NORM OF THE RESIDUAL

// 	//Defining res_h cd ..
// 	std::vector<std::vector<double > > res_h(grid_size_y, std::vector<double>(grid_size_x, 0.0));
// 
// 	computeRes(u,res_h, f);
// 	/*
// 	//set the res to be equal to 0 at the boundary
// 	for(std::size_t j=grid_size_quad-1; j<grid_size; ++j) {
// 	  res_h[grid_size_quad-1][j] = 0;
// 	}*/
// 
// 	//computing the discrete L2 norm
// 	double res_l2 = 0;
// 	for (std::size_t i=0;i<grid_size_y;++i){
// 		for (std::size_t j=0;j<grid_size_x;++j){
// 			res_l2 += 2* pow(res_h[i][j],2);
// 		}
// 	}
// 	res_l2 = res_l2/(grid_size*grid_size);
// 	res_l2 = sqrt(res_l2);
// 
// 	std::cout<<"the l2 norm of the residual is "<<res_l2<<"."<<std::endl;

	}
}


