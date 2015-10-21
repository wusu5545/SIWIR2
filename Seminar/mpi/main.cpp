#include "MGSolver.h"
#include <sys/time.h>
#include <fstream>

void Initialization(int l)
{
  MGSolver mg(l);
  Grid initl(l),initr(l);
  Grid f(l);
  real r,phi;
  
  ofstream initf;
  initf.open("init.dat");

  //Initialization
  mg.setGhostLayer(0,1);
  initl.resize(initl.ysize(),initl.xsize()+1);
  initr.resize(initr.ysize(),initr.xsize()+1);
  initl.fill(0.0);
  initr.fill(0.0);
  for (size_t i = 0; i< initl.ysize();++i){
	r = sqrt(1+pow((-1 + i*initl.hsize()),2));
	phi = acos(-1/r);
	initl(i,0) = sqrt(r)*sin(0.5*phi);
	phi = acos(1/r);
	initr(i,0) = sqrt(r)*sin(0.5*phi);
  }
  for (size_t j = 0; j< initl.xsize()-1;++j){
      r = sqrt(1+pow((-1 + j*initl.hsize()),2));
      phi = acos((-1+j*initl.hsize())/r);
      initl(0,j) = sqrt(r)*sin(0.5*phi);
  }
  for (size_t j = 0;j<initr.xsize()-1;++j){
      r = sqrt(1+pow((1 - j*initr.hsize()),2));
      phi = acos((1 - j*initr.hsize())/r);
      initr(0,j) = sqrt(r)*sin(0.5*phi);
  }
  
  for (size_t i=0;i<initl.ysize();++i){
    for (size_t j=0;j<initl.xsize()-1;++j)
	initf<<std::setw(10)<<(-1+j*initl.hsize())<<" "<<std::setw(10)<<(-1+i*initl.hsize())<<" "<<std::setw(10)<<initl(i,j)<<endl;
    for (size_t j=initr.xsize()-2;j>0;--j)
	initf<<std::setw(10)<<(1-(j-1)*initr.hsize())<<" "<<std::setw(10)<<(-1+i*initr.hsize())<<" "<<std::setw(10)<<initr(i,j-1)<<endl;
  }
  for (size_t i=initl.ysize()-1;i>0;--i){
    for (size_t j=0;j<initl.xsize()-1;++j)
	initf<<std::setw(10)<<(-1+j*initl.hsize())<<" "<<std::setw(10)<<(1-(i-1)*initl.hsize())<<" "<<std::setw(10)<<initl(i-1,j)<<endl;
    for (size_t j=initr.xsize()-2;j>0;--j)
	initf<<std::setw(10)<<(1-(j-1)*initr.hsize())<<" "<<std::setw(10)<<(1-(i-1)*initr.hsize())<<" "<<std::setw(10)<<initr(i-1,j-1)<<endl;    
  }

  initf.close();
  
  //mpi Initialization
  // Definition of the variables
  int size; //The total number of processes
  int rank; //The rank/number of this process
  // MPI initialization

  // Determining the number of CPUs and the rank for each CPU
  MPI_Comm_size( MPI_COMM_WORLD, &size );
  MPI_Comm_rank( MPI_COMM_WORLD, &rank );
  
  f.fill(0.0);
  
  std::cout<<"Your Alias: "<<"Wu&Lou"<<std::endl;
  struct timeval t0, t;
  gettimeofday(&t0, NULL);
  
  //pre smooth = 2, post smooth = 1 =>Vcycle(2,1)
  if( rank == 0 )
    for (size_t i = 0;i<8;++i)
      mg.solveFMG(initl,f,3,2,rank);
  if (rank == 1)
    for (size_t i = 0;i<8;++i)
      mg.solveFMG(initr,f,3,2,rank);

  gettimeofday(&t, NULL);
  std::cout << "Wall clock time of MG execution: " <<
  ((int64_t)(t.tv_sec - t0.tv_sec) * (int64_t)1000000 +
  (int64_t)t.tv_usec - (int64_t)t0.tv_usec) * 1e-3
  << " ms" << std::endl;  
  
  real send_buffer[initr.xsize()*initr.ysize()],recv_buffer[initr.xsize()*initr.ysize()];
  MPI_Request request;
  MPI_Status status;

  if(rank == 1){
    for (size_t i=1;i<initr.ysize()-1;++i)
      for (size_t j = 1;j<initr.xsize()-1;j++)
	send_buffer[i*initr.xsize()+j] = initr(i,j);
    MPI_Isend (send_buffer, initr.xsize()*initr.ysize(), MPI_DOUBLE, 0, 77, MPI_COMM_WORLD, &request);
    MPI_Wait (&request, &status);
  }
  if (rank == 0){
    MPI_Irecv (recv_buffer, initr.xsize()*initr.ysize(), MPI_DOUBLE, 1, 77, MPI_COMM_WORLD,&request);
    MPI_Wait (&request, &status);
    for (size_t i=1;i<initr.ysize()-1;++i)
      for (size_t j = 1;j<initr.xsize()-1;j++)
	initr(i,j) = recv_buffer[i*initr.xsize()+j];

    
  }
  
  // MPI finalizations
  MPI_Finalize();

  //output solution into "solution.txt"
  ofstream solution;
  solution.open("solution.dat");
  
  for (size_t i=0;i<initl.ysize();++i){
    for (size_t j=0;j<initl.xsize()-1;++j)
	solution<<std::setw(10)<<(-1+j*initl.hsize())<<" "<<std::setw(10)<<(-1+i*initl.hsize())<<" "<<std::setw(10)<<initl(i,j)<<endl;
    for (size_t j=initr.xsize()-2;j>0;--j)
	solution<<std::setw(10)<<(1-(j-1)*initr.hsize())<<" "<<std::setw(10)<<(-1+i*initr.hsize())<<" "<<std::setw(10)<<initr(i,j-1)<<endl;
  }
  for (size_t i=initl.ysize()-1;i>0;--i){
    for (size_t j=0;j<initl.xsize()-1;++j)
	solution<<std::setw(10)<<(-1+j*initl.hsize())<<" "<<std::setw(10)<<(1-(i-1)*initl.hsize())<<" "<<std::setw(10)<<initl(i-1,j)<<endl;
    for (size_t j=initr.xsize()-2;j>0;--j)
	solution<<std::setw(10)<<(1-(j-1)*initr.hsize())<<" "<<std::setw(10)<<(1-(i-1)*initr.hsize())<<" "<<std::setw(10)<<initr(i-1,j-1)<<endl;    
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
  MPI_Init( &argc, &argv);
  Initialization(atoi(argv[1]));
  return 0;
}