#include <math.h>
#include <iostream>
#include <fstream>
#include "FUNCTIONS.h"
#include <vector>
#include <iomanip>

void fmg (std::vector<std::vector<double > > &u, std::vector<std::vector<double > > &f,int v1, int v2, int gamma, int miu,
	 int l, std::size_t grid_size)
{

      std::size_t grid_size_y = u.size();
      std::size_t grid_size_x = u[0].size();

      if(l==2) {
	rbgs(u, 10, f);
      }
  
      else
    {
      std::size_t coarser_grid_size_x = (grid_size_x-1)/2+1;
      std::size_t coarser_grid_size_y = (grid_size_y-1)/2+1;

      
      std::vector<std::vector<double > > res_h(grid_size_y, std::vector<double>(grid_size_x, 0.0));
      
      computeRes(u, res_h, f);
           
      std::vector<std::vector<double > > f_2h(coarser_grid_size_y, std::vector<double>(coarser_grid_size_x, 0.0));
      std::vector<std::vector<double > > c_2h(coarser_grid_size_y, std::vector<double>(coarser_grid_size_x, 0.0));

      restriction(res_h, f_2h);
      
      fmg(c_2h, f_2h, v1, v2, gamma, miu, l-1, coarser_grid_size_x);


      correction(u, c_2h);

      for (int z =0; z<miu; ++z) {
	mgm(u,f, v1, v2, gamma,l, grid_size);
      }
  }

}