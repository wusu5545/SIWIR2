#include <math.h>
#include <iostream>
#include <fstream>
#include "FUNCTIONS.h"
#include <vector>
#include <iomanip>

void mgm(std::vector<std::vector<double > > &u, std::vector<std::vector<double > > &f,int v1, int v2, int gamma, int l, std::size_t grid_size) {

  std::size_t grid_size_y = u.size();
  std::size_t grid_size_x = u[0].size();

  //PERFORM V1 PRE-SMOOTHERNING STEPS
  rbgs(u, v1, f);

/*
    solution<<"the velocity distribution after the first relaxation:"<< std::endl;

  for (int i=0; i<grid_size_y; ++i){
    for(int j=0; j<grid_size_x; ++j){
      solution<<std::setw(10)<<u[i][j]<<" ";
    }
    solution<<std::endl;
  }
*/

  //COMPUTE RES

  std::vector<std::vector<double > > res_h(grid_size_y, std::vector<double>(grid_size_x, 0.0));

  computeRes(u, res_h, f);
/*
  solution<<"the residual after the first relaxation:"<< std::endl;
  for (std::size_t i=0; i<grid_size_y; ++i){
    for(std::size_t j=0; j<grid_size_x; ++j){
      solution<<std::setw(10)<<res_h[i][j]<<" ";
    }
    solution<<std::endl;
  }*/

  //RESTRICTION OF RESIDUAL
  std::size_t coarser_grid_size_x = (grid_size_x-1)/2+1;
  std::size_t coarser_grid_size_y = (grid_size_y-1)/2+1;

  std::vector<std::vector<double > > f_2h(coarser_grid_size_y, std::vector<double>(coarser_grid_size_x, 0.0));

  restriction(res_h, f_2h);
/*
   solution<<"the f_2h after the restriction:"<< std::endl;
  for (std::size_t i=0; i<coarser_grid_size_y; ++i){
    for(std::size_t j=0; j<coarser_grid_size_x; ++j){
      solution<<std::setw(10)<<f_2h[i][j]<<" ";
    }
    solution<<std::endl;
  }*/
    //SET UP c_2h

 std::vector<std::vector<double > > c_2h(coarser_grid_size_y, std::vector<double>(coarser_grid_size_x, 0.0));

  //SOLVE IF LEVEL

  if (l==3){
    rbgs(c_2h, 7, f_2h);
/*
        solution<<"the correction after solving is:"<< std::endl;
  for (int i=0; i<coarser_grid_size_y; ++i){
    for(int j=0; j<coarser_grid_size_x; ++j){
      solution<<std::setw(10)<<c_2h[i][j]<<" ";
    }
    solution<<std::endl;
  }
*/
   }

  else{
    for(int k=0; k<gamma; ++k)
    {
      mgm(c_2h,f_2h, v1, v2, gamma,l-1, grid_size);
    }
  }


  //CORRECT
  correction(u, c_2h);
/*
         solution<<"the velocity distribution after correction:"<< std::endl;
  for (int i=0; i<grid_size_y; ++i){
    for(int j=0; j<grid_size_x; ++j){
      solution<<std::setw(10)<<u[i][j]<<" ";
    }
    solution<<std::endl;
  }
*/

  //PERFORM v2 POST SMOOTHERNING STEPS
  rbgs(u, v2, f);


  //computeRes(u, res_h,f);
/*
       solution<<"the residual after the post relaxation:"<< std::endl;
  for (int i=0; i<grid_size_y; ++i){
    for(int j=0; j<grid_size_x; ++j){
      solution<<std::setw(10)<<res_h[i][j]<<" ";
    }
    solution<<std::endl;
  }*/

}



