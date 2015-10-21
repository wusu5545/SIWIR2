#include "functions.h"
#include <sys/time.h>

#include <boost/shared_ptr.hpp>
#include <boost/make_shared.hpp>

using namespace lbm;

void Solver(FileReader &reader)
{
  //Converting the 'timesteps' parameter to size_t value
  const string geometry(
    reader.getFileName<string>("geometry"));
  
  bool CASE;
  size_t sizex=0,sizey=0;
  size_t vtk_count = 0;
  Flags  fluid(sizex,sizey);
  
  if (geometry == ""){
    CASE = 0;
    sizex = reader.getParameter<size_t>("sizex");
    sizey = reader.getParameter<size_t>("sizey");
  }else
  {
    CASE = 1;
    fstream input;
    input.open(geometry,fstream::in);
    string inputLine = "";
    
    // First line : version
    getline(input,inputLine);
    if(inputLine.compare("P2") != 0) 
      cerr << "Version error" << endl;
    else 
      cout << "Version : " << inputLine << endl;
    
    // Second line : comment
    getline(input,inputLine);
    cout << "Comment : " << inputLine << endl;
    
    // Third line : size
    input >> sizex >> sizey;
    cout << sizex << " columns and " << sizey << " rows" << endl;
    fluid.resize(sizex+2,sizey+2);

    // Following lines : data
    size_t pixel;
    // biggest number
    input >> pixel;
    
    for(size_t j = 1; j < sizey+1; ++j)
      for (size_t i = 1; i < sizex+1; ++i){ 
	input >> pixel;
	if (pixel == 255) 
	  fluid(i,sizey+1-j) = 1;
	else
	  fluid(i,sizey+1-j) = 0;
      }

    // Now print the array to see the result
    for(size_t j = 1; j < sizey+1; ++j) {
      for(size_t i = 1; i < sizex+1; ++i) {
	cout << fluid(i,sizey+1-j) << " ";
      }
      cout << endl;
    }
    input.close();
  }
  const size_t timesteps(
    reader.getParameter<size_t>("timesteps"));
  const double omega(
    reader.getParameter<double>("omega"));
  const size_t vtk_step(
    reader.getParameter<size_t>("vtk_step"));
  string vtk_file(
    reader.getFileName<string>("vtk_file"));
  
  vtk_file = string(vtk_file.begin(),vtk_file.end()-4);

  //Checking the value
  if (timesteps == 0 || timesteps > 1000000){
    cerr<< " Invalid 'timesteps' parameter!"<<endl;
    exit(0);
  }

  size_t Num_Fluid_Cells = sizex * sizey;
  
  sizex += 2;
  sizey += 2;
  
  //specifying Grids
  PDF_Field f_src(sizex,sizey);
  PDF_Field f_dst(sizex,sizey);
  PDF_Field f_eqm(sizex,sizey);
  
  V_Field   v(sizex,sizey);
  D_Field   rho(sizex,sizey);
  
  Flags	    flags(sizex,sizey);
  
  V_Field   c_alpha(9,1);//Cell
  D_Field   t_alpha(9,1);
  
  //Initialization
  c_alpha(0,0,0) =  0.0; c_alpha(0,0,1) =  0.0;
  c_alpha(1,0,0) =  0.0; c_alpha(1,0,1) =  1.0;
  c_alpha(2,0,0) =  0.0; c_alpha(2,0,1) = -1.0;
  c_alpha(3,0,0) = -1.0; c_alpha(3,0,1) =  0.0;
  c_alpha(4,0,0) =  1.0; c_alpha(4,0,1) =  0.0;
  c_alpha(5,0,0) = -1.0; c_alpha(5,0,1) =  1.0;
  c_alpha(6,0,0) =  1.0; c_alpha(6,0,1) =  1.0;
  c_alpha(7,0,0) = -1.0; c_alpha(7,0,1) = -1.0;
  c_alpha(8,0,0) =  1.0; c_alpha(8,0,1) = -1.0;
/* 5   1   6
    \  |  /
     \ | /
      \|/
  3----0----4  
      /|\
     / | \
    /  |  \
   7   2   8	*/

  t_alpha(0,0) = 4.0/9.0;
  t_alpha(1,0) = 1.0/9.0;
  t_alpha(2,0) = 1.0/9.0;
  t_alpha(3,0) = 1.0/9.0;
  t_alpha(4,0) = 1.0/9.0;
  t_alpha(5,0) = 1.0/36.0;
  t_alpha(6,0) = 1.0/36.0;
  t_alpha(7,0) = 1.0/36.0;
  t_alpha(8,0) = 1.0/36.0;
  
  //flags
  for (size_t j = 1;j < sizey - 1;++j)
    for (size_t i = 1;i <sizex - 1;++i)
      flags(i,j) = 1;
  if (CASE == 0){
    fluid.resize(sizex,sizey);
    for (size_t j = 1;j < sizey - 1;++j)
      for (size_t i = 1;i <sizex - 1;++i)
	fluid(i,j) = 1;
  }

  //density
  for (size_t j = 1;j < sizey - 1;++j)
    for (size_t i = 1;i <sizex - 1;++i)
      if (fluid(i,j) == 1)
	rho(i,j) = 1.0;
  //velocity
  for (size_t i = 0; i<sizex; ++i)
    v(i,sizey-1,0) = 0.08;
  
  //PDF
  for (size_t j = 1; j<sizey - 1;++j)
    for (size_t i = 1;i<sizex - 1;++i)
      if (fluid(i,j) == 1)
	for (size_t k = 0;k<9;++k)
	  f_src(i,j,k) = t_alpha(k,0);
      
    
  struct timeval t0, t;
  gettimeofday(&t0, NULL);
  
  for (size_t t = 0;t<=timesteps;++t){
    //Treat Boundary
    for (size_t j = 0;j<sizey;++j){
      if (j!=0){
	f_src(0,j,8) = f_src(1,j-1,5);
	f_src(sizex-1,j,7) = f_src(sizex-2,j-1,6); 
      }
      f_src(0,j,4) = f_src(1,j,3);
      f_src(sizex-1,j,3) = f_src(sizex-2,j,4);
      if (j+1!=sizey){
	f_src(0,j,6) = f_src(1,j+1,7);
	f_src(sizex - 1,j,5) = f_src(sizex-2,j+1,8);
      }
    }
    for (size_t i=0;i<sizex;++i){
      if(i!=0){
	f_src(i,0,5) = f_src(i-1,1,8);
	f_src(i,sizey-1,7) = f_src(i-1,sizey-2,6) - 
		  2*t_alpha(6,0)*3*(c_alpha(6,0,0)*v(i,sizey-1,0)+c_alpha(6,0,1)*v(i,sizey-1,1));
      }
      f_src(i,0,1) = f_src(i,1,2);
      f_src(i,sizey-1,2) = f_src(i,sizey-2,1) - 
		    2*t_alpha(1,0)*3*(c_alpha(1,0,0)*v(i,sizey-1,0)+c_alpha(1,0,1)*v(i,sizey-1,1));
      if(i+1!=sizex){
	f_src(i,0,6) = f_src(i+1,1,7);
	f_src(i,sizey-1,8) = f_src(i+1,sizey-2,5) - 
		  2*t_alpha(5,0)*3*(c_alpha(5,0,0)*v(i,sizey-1,0)+c_alpha(5,0,1)*v(i,sizey-1,1)); 
      }
    }
    
    if (CASE == 1)
      for (size_t j=1;j<sizey-1;++j)
	for (size_t i=1;i<sizex-1;++i)
	  if (fluid(i,j) == 0){
	    f_src(i,j,0) = f_src(i  ,j  ,0);
	    f_src(i,j,1) = f_src(i  ,j+1,2);
	    f_src(i,j,2) = f_src(i  ,j-1,1);
	    f_src(i,j,3) = f_src(i-1,j  ,4);
	    f_src(i,j,4) = f_src(i+1,j  ,3);
	    f_src(i,j,5) = f_src(i-1,j+1,8);
	    f_src(i,j,6) = f_src(i+1,j+1,7);
	    f_src(i,j,7) = f_src(i-1,j-1,6);
	    f_src(i,j,8) = f_src(i+1,j-1,5);
	  }
	    
    
    //Streaming step (PULL)
    for(size_t j = 1;j<sizey-1;++j)
      for (size_t i = 1;i<sizex-1;++i)
	if (fluid(i,j) == 1)
	{
	  f_dst(i,j,0) = f_src(i  ,j  ,0);
	  f_dst(i,j,1) = f_src(i  ,j-1,1);
	  f_dst(i,j,2) = f_src(i  ,j+1,2);
	  f_dst(i,j,3) = f_src(i+1,j  ,3);
	  f_dst(i,j,4) = f_src(i-1,j  ,4);
	  f_dst(i,j,5) = f_src(i+1,j-1,5);
	  f_dst(i,j,6) = f_src(i-1,j-1,6);
	  f_dst(i,j,7) = f_src(i+1,j+1,7);
	  f_dst(i,j,8) = f_src(i-1,j+1,8);
	}
    
    //Collide
    //Calculation of Macroscopic density and velocity
    for(size_t j=1;j<sizey-1;++j)
      for(size_t i=1;i<sizex-1;++i){
	rho(i,j) = 0;
	for (size_t k = 0;k<9;++k)
	  rho(i,j) += f_dst(i,j,k);
      }
    
    for (size_t j=1;j<sizey-1;++j)
      for (size_t i=1;i<sizex-1;++i){
	v(i,j,0) = 0;
	v(i,j,1) = 0;
	if (fluid(i,j) == 1)
	  for (size_t k=0;k<9;++k){
	    v(i,j,0) += f_dst(i,j,k)*c_alpha(k,0,0);
	    v(i,j,1) += f_dst(i,j,k)*c_alpha(k,0,1);
	  }
//  	v(i,j,0) = v(i,j,0)/rho(i,j);
//  	v(i,j,1) = v(i,j,1)/rho(i,j);
      }
    
    //Calculation of Equilibrium Distribution Function
    for (size_t j=1;j<sizey-1;++j)
      for (size_t i=1;i<sizex-1;++i)
	if (fluid(i,j) == 1)
	  for (size_t k=0;k<9;++k)
	    f_eqm(i,j,k) = t_alpha(k,0)*(rho(i,j)+3.0*(c_alpha(k,0,0)*v(i,j,0)+c_alpha(k,0,1)*v(i,j,1))+
			  9.0/2.0*((c_alpha(k,0,0)*v(i,j,0)+c_alpha(k,0,1)*v(i,j,1))*(c_alpha(k,0,0)*v(i,j,0)+c_alpha(k,0,1)*v(i,j,1)))-
			  3.0/2.0*(v(i,j,0)*v(i,j,0)+v(i,j,1)*v(i,j,1)));
    
    //Collistion step
    for (size_t j=1;j<sizey-1;++j)
      for (size_t i=1;i<sizex-1;++i)
	if (fluid(i,j) == 1)
	  for (size_t k=0;k<9;++k)
	    f_dst(i,j,k) = f_dst(i,j,k) - omega * (f_dst(i,j,k) - f_eqm(i,j,k));
	
    f_src.swap(f_dst);
    f_dst.reset();
    
    vtk_output(vtk_step,t,vtk_file,vtk_count,sizex,sizey,flags,rho,v);
  }
  
   //FINAL MEASURE OF PERFORMANCE
 
  gettimeofday(&t, NULL);
  double total_time = ((int64_t)(t.tv_sec - t0.tv_sec) * (int64_t) 1000000 + (int64_t)t.tv_usec - (int64_t)t0.tv_usec) * 1e-6;
  double MLUps = (Num_Fluid_Cells*timesteps)/total_time;

  cout<<"Mega Lattice site Updates per second is "<<MLUps<<endl;
  
}