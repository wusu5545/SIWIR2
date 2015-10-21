#include <math.h>
#include  <iostream>
#include <fstream>
#include "FUNCTIONS.h"
#include <vector>
#include <iomanip>
#include <stdlib.h>
#include <sys/time.h>
#include <omp.h>

int main(int argc, char *argv[])
{
  std::cout<<"Your Alias: "<< "Ang_Wang"<<std::endl;
  
  
  //READING IN FROM TERMINAL

  int l;
  if(argc!=1)
    l = atoi(argv[1]);
  else
    l = 11;

   //FOR OUTPUT
    std::ofstream init;
    init.open("init.dat");

    std::ofstream solution;
    solution.open("solution.dat");

   //DEFINITION OF PARAMETERS
    std::size_t grid_size = pow(2,l)+1;

    double h = 2.0/((double) (grid_size-1));

    std::size_t grid_size_x = grid_size;
    std::size_t grid_size_y = (grid_size - 1)/2 + 1;

    //STORING AND PRINTING THE ENTIRE VELOCITY DIST
    std::vector<std::vector<double > > u(grid_size_y, std::vector<double>(grid_size_x));

    init_and_boundary(u);

    std::vector<std::vector<double > > f(grid_size_y, std::vector<double>(grid_size_x));

    for (std::size_t i=0;i<grid_size_y;++i){
	for (std::size_t j=0;j<grid_size_x;++j){
		init<<std::setw(10)<<(j*h-1)<<" "<<std::setw(10)<<(1-i*h)<<" "<<std::setw(10)<<u[i][j]<<std::endl;
	}
    }
    for (std::size_t i=0; i<grid_size_y-1 ;++i){
	for (std::size_t j=0;j<grid_size_x;++j){
		init<<std::setw(10)<<(j*h-1)<<" "<<std::setw(10)<<(-1+i*h)<<" "<<std::setw(10)<<u[i][j]<<std::endl;
	}
    }
    init.close();
    
    
  //TO MEASURE THE WALL CLOCK TIME
  
  struct timeval t0, t;

  gettimeofday(&t0, NULL);

 solveMG(l, u, f);
 
  //MEASURE WALL CLOCK TIME THE PROGRAMME TAKES
 gettimeofday(&t, NULL);

  std::cout<<"Wall clock time of MG execution: "<< ((int64_t)(t.tv_sec - t0.tv_sec) * (int64_t) 1000000 + (int64_t)t.tv_usec - (int64_t)t0.tv_usec) * 1e-3 << " ms"<<std::endl;

  // WRITING THE COMPUTED SOLUTION TO A FILE NAMED SOLUTION.TXT

     for (std::size_t i=0;i<grid_size_y;++i){
	for (std::size_t j=0;j<grid_size_x;++j){
		solution<<std::setw(10)<<(j*h-1)<<" "<<std::setw(10)<<(1-i*h)<<" "<<std::setw(10)<<u[i][j]<<std::endl;
	}
  }
    for (std::size_t i=0; i<grid_size_y-1 ;++i){
	for (std::size_t j=0;j<grid_size_x;++j){
		solution<<std::setw(10)<<(j*h-1)<<" "<<std::setw(10)<<(-1+i*h)<<" "<<std::setw(10)<<u[i][j]<<std::endl;
	}
  }   
   solution.close();
}



