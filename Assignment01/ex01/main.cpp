#include "MGSolver.h"
#include <sys/time.h>
#include <fstream>

#define CASE 0 //0 stands for the main assign,1 stands for the bonus

void Initialization(int l, int n)
{
  MGSolver mg(l);
  Grid init(l);
  Grid f(l);
  if (CASE == 0){
    real tol = pow(10.0,-10);
    for (size_t i = 0; i< init.ysize();i++){
	  init(i,0) = sin(M_PI*(i)*init.hsize())*sinh(M_PI*0);
	  init(i,init.xsize()-1) = sin(M_PI*(i)*init.hsize())*sinh(M_PI*(1));
	  if(init(i,0) < tol){init(i,0) = 0;}
	  if(init(i,init.xsize()-1) < tol){init(i,init.xsize()-1) = 0;}
    }
    for (size_t j = 0; j< init.xsize();j++){
	init(0,j) = sin(M_PI*0)*sinh(M_PI*(j)*init.hsize());
	init(init.ysize()-1,j) = sin(M_PI*(1))*sinh(M_PI*(j)*init.hsize());
	if(init(0,j) < tol){init(0,j) = 0;}
	if(init(init.ysize()-1,j) < tol){init(init.ysize()-1,j) = 0;}
    }
  }
  else{
    init.resize(init.ysize(),init.xsize()+2);
    f.resize(f.ysize(),f.xsize()+2);
    init.fill(0.0);
    mg.setGhostLayer(0,2);
    for (size_t i = 1;i<init.ysize()-1;++i){
      init(0,i+1) = i*init.hsize()*(1-i*init.hsize());
      init(init.ysize()-1,i+1) = i*init.hsize()*(1-i*init.hsize());
    }
    f.fill(2.0);
  }
  
  struct timeval start,end;
  
  gettimeofday(&start,NULL);
  
  //pre smooth = 2, post smooth = 1 =>Vcycle(2,1)
  mg.solve(init,f,n,2,1);
  
  gettimeofday(&end,NULL);
  
  real seconds,useconds, total_time;
  
  seconds = end.tv_sec - start.tv_sec;
  useconds = end.tv_usec - start.tv_usec;
  
  total_time = (seconds*1000.0 + useconds/1000.0);
  
  cout<<"The total wall clock time of "<<n<<" V-cycles is "<<total_time<<"ms"<<endl;

  //output solution into "solution.txt" and "error.txt"
  ofstream solution;
  solution.open("solution.txt");
  
  solution << "#" <<setw(12)<<"x"<<setw(12)<<"y"<<setw(12)<<"u"<<endl;
  if (CASE == 0)
    for (size_t i=0;i<init.ysize();++i)
      for (size_t j=0;j<init.xsize();++j)
	solution<<setw(12)<<i*init.hsize()<<" "<<setw(12)<<j*init.hsize()<<setw(12)<<init(i,j)<<endl;
  else
    for (size_t i=0;i<init.ysize();++i)
      for (size_t j=1;j<init.xsize()-1;j++)
      {
	//init(i,j) = (1-i*init.hsize())*i*init.hsize();
	solution<<setw(12)<<i*init.hsize()<<" "<<setw(12)<<j*init.hsize()<<setw(12)<<init(i,j)<<endl;
      }
  solution.close();

  ofstream error;
  error.open("error.txt",std::fstream::app);

  //measure the error in L2 norm
  real err = 0;
  if (CASE == 0)
    for (size_t i = 0;i<init.ysize();i++){
      for (size_t j = 0;j<init.xsize();j++){
	err += pow(init(i,j)-sin(M_PI*i*init.hsize())*sinh(M_PI*j*init.hsize()),2);
      }
    }
  else
    for (size_t i = 0;i<init.ysize();i++){
      for (size_t j = 1;j<init.xsize()-1;j++){
	err += pow(init(i,j)-(1-i*init.hsize())*i*init.hsize(),2);
      }
    }
  
  err = sqrt(err/(init.xsize()*init.ysize()));

  error<<CASE<<setw(12)<<n<<setw(12)<<init.hsize()<<" "<<setw(12)<<err<<endl;
  
  error.close();
  
}

int main(int argc, char** argv)
{
  if (argc !=3)
  {
    cout<< "Args: levels and cycles" <<endl;
    exit(0);
  }
  Initialization(atoi(argv[1]),atoi(argv[2]));
  return 0;
}
