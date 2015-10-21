#include "MGSolver.h"
#include <sys/time.h>
#include <fstream>


void Initialization(int l)
{
  MGSolver mg(l);
  Grid init(l);
  Grid f(l);
  real r,phi;
  
  init.fill(0.0);
  ofstream initf;
  initf.open("init.dat");

  //Initialization
  init.fill(0.0);
  for (size_t i = 0; i< init.ysize();i++){
	r = sqrt(1+pow((-1 + i*init.hsize()),2));
	phi = acos(-1/r);
	init(i,0) = sqrt(r)*sin(0.5*phi);
	phi = acos(1/r);
	init(i,init.xsize()-1) = sqrt(r)*sin(0.5*phi);
  }
  for (size_t j = 0; j< init.xsize();j++){
      r = sqrt(1+pow((-1 + j*init.hsize()),2));
      phi = acos((-1+j*init.hsize())/r);
      init(0,j) = sqrt(r)*sin(0.5*phi);
  }
  
  for (size_t i=0;i<init.ysize();++i)
    for (size_t j=0;j<init.xsize();++j){
	initf<<std::setw(10)<<(-1+j*init.hsize())<<" "<<std::setw(10)<<(-1+i*init.hsize())<<" "<<std::setw(10)<<init(i,j)<<endl;
    }
  for (size_t i=init.ysize()-1;i>0;--i)
    for (size_t j=0;j<init.xsize();++j){
	initf<<std::setw(10)<<(-1+j*init.hsize())<<" "<<std::setw(10)<<(1-(i-1)*init.hsize())<<" "<<std::setw(10)<<init(i-1,j)<<endl;
    }

  initf.close();
  f.fill(0.0);
  
  std::cout<<"Your Alias: "<<"Wu&Lou"<<std::endl;
  struct timeval t0, t;
  gettimeofday(&t0, NULL);
  
  //pre smooth = 2, post smooth = 1 =>Vcycle(2,1)
  for (size_t i = 0;i<1;++i){
    mg.solveFMG(init,f,1,1);
  }
  //mg.solve(init,f,1,3,2);
  
  gettimeofday(&t, NULL);
  std::cout << "Wall clock time of MG execution: " <<
  ((int64_t)(t.tv_sec - t0.tv_sec) * (int64_t)1000000 +
  (int64_t)t.tv_usec - (int64_t)t0.tv_usec) * 1e-3
  << " ms" << std::endl;  

  //output solution into "solution.txt" and "error.txt"
  ofstream solution;
  solution.open("solution.dat");
  
  for (size_t i=0;i<init.ysize();++i)
    for (size_t j=0;j<init.xsize();++j){
	solution<<std::setw(10)<<(-1+j*init.hsize())<<" "<<std::setw(10)<<(-1+i*init.hsize())<<" "<<std::setw(10)<<init(i,j)<<endl;
    }
  for (size_t i=init.ysize()-1;i>0;--i)
    for (size_t j=0;j<init.xsize();++j){
	solution<<std::setw(10)<<(-1+j*init.hsize())<<" "<<std::setw(10)<<(1-(i-1)*init.hsize())<<" "<<std::setw(10)<<init(i-1,j)<<endl;
    }

  solution.close();
  
}

int main(int argc, char** argv)
{
  if (argc !=2)
  {
    cout<< "Args: levels" <<endl;
    exit(0);
  }
  Initialization(atoi(argv[1]));
  return 0;
}