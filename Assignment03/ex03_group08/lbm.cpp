#include "FUNCTIONS.h"
#include "GRID.h"
#include <fstream>
#include <iostream>
#include <cmath>
#include <iomanip>
#include <string>
#include<sstream>
#include <fstream>
#include <sys/time.h>

using namespace lbm;

int main(int argc, char *argv[]) {
  
  struct timeval t0, t;

  gettimeofday(&t0, NULL);
  
  if (argc < 2)
  {
    std::cout << "No input file specified" << std::endl;
    return 0;
  }
  
  int xsize;
  int ysize;
  uint timesteps;
  double omega;
  int vtk_step;
  std::string vtk_file_long;
  int vtk_out_count = 0;
  
  
  //READING INPUT
  std::string temp;
  std::string par_file_name (argv[1]);
  std::ifstream par_file(par_file_name.c_str());
  
  while(!par_file.eof())
  {
    par_file >> temp;

		if(temp.compare("sizex") == 0)
			par_file >> xsize;
		else if(temp.compare("sizey") == 0)
			par_file >> ysize;
		else if(temp.compare("timesteps") == 0)
			par_file >> timesteps;
		else if(temp.compare("omega") == 0)
			par_file >> omega;
		else if(temp.compare("vtk_file") == 0)
			par_file >> vtk_file_long;
		else if(temp.compare("vtk_step") == 0)
			par_file >> vtk_step;
  }

  xsize +=2;
  ysize +=2; 
  std::string vtk_file (vtk_file_long.begin(), vtk_file_long.end()-4); 
  
  //FOR PERFORMANCE
  int num_fluid_cells = (xsize-2)*(ysize-2);
  double MLUps = 0;


  //SPECIFYING ALL THE GRIDS NEEDED
  PDF_Field f_source(xsize, ysize);
  PDF_Field f_dest(xsize, ysize);
  PDF_Field f_eqm(xsize, ysize);
  V_Field v(xsize, ysize);
  D_Field d(xsize, ysize);
  FLAGS flags(xsize, ysize);
  Grid<double, 2> c_alpha(9,1);

  c_alpha(0,0,0) = 0.0; c_alpha(0,0,1) = 0.0;
  c_alpha(1,0,0) = 0.0; c_alpha(1,0,1) = 1.0;
  c_alpha(2,0,0) = 0.0; c_alpha(2,0,1) = -1.0;
  c_alpha(3,0,0) = -1.0; c_alpha(3,0,1) = 0.0;
  c_alpha(4,0,0) = 1.0; c_alpha(4,0,1) = 0.0;
  c_alpha(5,0,0) = -1.0; c_alpha(5,0,1) = 1.0;
  c_alpha(6,0,0) = 1.0; c_alpha(6,0,1) = 1.0;
  c_alpha(7,0,0) = -1.0; c_alpha(7,0,1) = -1.0;
  c_alpha(8,0,0) = 1.0; c_alpha(8,0,1) = -1.0;

  Grid<double, 1> t_alpha(9,1);
  t_alpha(0,0) = 4.0/9.0;
  t_alpha(1,0) = 1.0/9.0;
  t_alpha(2,0) = 1.0/9.0;
  t_alpha(3,0) = 1.0/9.0;
  t_alpha(4,0) = 1.0/9.0;  
  t_alpha(5,0) = 1.0/36.0;
  t_alpha(6,0) = 1.0/36.0;
  t_alpha(7,0) = 1.0/36.0;
  t_alpha(8,0) = 1.0/36.0;

  //INITIALIZATION
  //flags
  for(int j=1; j<ysize-1; ++j){
  for(int i=1; i<xsize-1; ++i){
    flags(i, j) = 1;
  }
  }

  //density
  for(int j=1; j<ysize-1; ++j){
    for(int i=1; i<xsize-1; ++i){
      d(i,j) = 1.0;
    }
  }
  
  //velocity
  for(int i=0; i<xsize; ++i) {
    v(i, ysize-1, 0) = 0.08;
  }

  
  //pdf
  for(int j=1; j<ysize-1; ++j){
    for(int i=1; i<xsize-1; ++i){
      for(int k=0; k<9; ++k){
	f_source(i,j,k) = t_alpha(k,0);
      }
    }
  }  
  
for(uint t = 0; t<timesteps;++t){

  //TREAT BOUNDARY
  for(int j=0; j<ysize; ++j){
    if(j+1!=ysize){
    f_source(0, j, 6) = f_source(1,j+1,7);}
    f_source(0, j, 4) = f_source(1,j,3);
    if(j-1!=-1){
    f_source(0, j, 8) = f_source(1,j-1,5); }
    
    if(j+1!=ysize){
    f_source(xsize-1, j, 5) = f_source(xsize-2, j+1,8); }
    f_source(xsize-1, j, 3) = f_source(xsize-2, j,4);
    if(j-1!=-1){
    f_source(xsize-1, j, 7) = f_source(xsize-2, j-1,6);}
  }
  for(int i=0; i<xsize; ++i){
    if(i-1!=-1){
    f_source(i, 0, 5) = f_source(i-1,1,8);}
    f_source(i, 0, 1) = f_source(i,1,2);
    if(i+1!=xsize){
    f_source(i, 0, 6) = f_source(i+1,1,7);} 
    
    if(i-1!=-1){
    f_source(i, ysize-1, 7) = f_source(i-1, ysize-2, 6) - 2*t_alpha(6,0)*3*(c_alpha(6,0,0)*v(i, ysize-1,0)+
			     c_alpha(6,0,1)*v(i, ysize-1,1));}
    f_source(i, ysize-1, 2) = f_source(i, ysize-2, 1) - 2*t_alpha(1,0)*3*(c_alpha(1,0,0)*v(i, ysize-1,0)+
 			     c_alpha(1,0,1)*v(i, ysize-1,1));
    if(i+1!=xsize){
    f_source(i, ysize-1, 8) = f_source(i+1, ysize-2, 5) - 2*t_alpha(5,0)*3*(c_alpha(5,0,0)*v(i, ysize-1,0)+
 			     c_alpha(5,0,1)*v(i, ysize-1,1));} 
  }  
  
    
   //STREAMING STEP (pull)
     for(int j=1; j<ysize-1; ++j){
    for(int i=1; i<xsize-1; ++i){
      f_dest(i,j,0) = f_source(i,j,0);
      f_dest(i,j,1) = f_source(i,j-1,1);
      f_dest(i,j,2) = f_source(i,j+1,2);
      f_dest(i,j,3) = f_source(i+1,j,3);
      f_dest(i,j,4) = f_source(i-1,j,4);
      f_dest(i,j,5) = f_source(i+1,j-1,5);
      f_dest(i,j,6) = f_source(i-1,j-1,6);
      f_dest(i,j,7) = f_source(i+1,j+1,7);
      f_dest(i,j,8) = f_source(i-1,j+1,8);
      }
    }

    // COLLIDE
    //CALCULATION OF MACROSCOPIC DENSITY AND VELOCITY

   for(int j=1; j<ysize-1; ++j){
    for(int i=1; i<xsize-1; ++i){
      d(i,j) = 0;
      for(int k=0; k<9; ++k) {
	d(i, j) += f_dest(i,j,k);
      }
    }
   }

   for(int j=1; j<ysize-1; ++j){
    for(int i=1; i<xsize-1; ++i){
      v(i,j,0) = 0;
      v(i,j,1) = 0;
      for(int k=0; k<9; ++k) {
	v(i, j, 0) += f_dest(i,j,k)*c_alpha(k,0,0);
	v(i, j, 1) += f_dest(i,j,k)*c_alpha(k,0,1);
      }
      //v(i,j,0) = v(i,j,0)/d(i,j);
      //v(i,j,1) = v(i,j,1)/d(i,j);
    }
   }

  //CALCULATION OF EQUILIBRIUM DISTRIBUTION FUNCTION
  for(int j=1; j<ysize-1; ++j){
    for(int i=1; i<xsize-1; ++i){
      for(int k=0; k<9; ++k){
	f_eqm(i,j,k) = t_alpha(k,0)*(d(i,j)+3.0*(c_alpha(k,0,0)*v(i,j,0) + c_alpha(k,0,1)*v(i,j,1))
	+ 9.0/2.0*pow(c_alpha(k,0,0)*v(i,j,0) + c_alpha(k,0,1)*v(i,j,1),2) - 3.0/2.0*(v(i,j,0)*v(i,j,0)
	+v(i,j,1)*v(i,j,1)));
      }
    }
  }

  //COLLISION STEP
   for(int j=1; j<ysize-1; ++j){
    for(int i=1; i<xsize-1; ++i){
      for(int k=0; k<9; ++k){
	f_dest(i,j,k) = f_dest(i,j,k) - omega*(f_dest(i,j,k) - f_eqm(i,j,k));
      }
    }
  }

  f_source.swap(f_dest);
  f_dest.reset();

  vtk_output(vtk_step, t, vtk_file, vtk_out_count,xsize, ysize, flags, d, v);
 }
 
 //FINAL MEASURE OF PERFORMANCE
 
  gettimeofday(&t, NULL);
  double total_time = ((int64_t)(t.tv_sec - t0.tv_sec) * (int64_t) 1000000 + (int64_t)t.tv_usec - (int64_t)t0.tv_usec) * 1e-6;
  MLUps = (num_fluid_cells*timesteps)/total_time;

  std::cout<<"Mega Lattice site Updates per second is "<<MLUps<<std::endl;

  return 0;
}