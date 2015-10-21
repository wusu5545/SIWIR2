#include <vector>
#include <cmath>
#include <iostream>
#include <omp.h>



void init_and_boundary(std::vector<std::vector <double> > & u){

  std::size_t grid_size_y = u.size();
  std::size_t grid_size_x = u[0].size();
   double h = 2.0/(double)(grid_size_x-1);
#pragma omp parallel for
    for(std::size_t j=0; j<grid_size_x; ++j){
	double r = sqrt((-1 + j*h)*(-1 + j*h) + 1);
	double	phi = acos((-1 + j*h) / r);
		u[0][j] = sqrt(r)*sin(phi / 2.0);
		//u[grid_size_y-1][j]= sqrt(r)*sin(phi / 2.0);
    }
#pragma omp parallel for
     for(std::size_t i=0; i<grid_size_y; ++i){
       double r = sqrt(1 + (1-i*h)*(1-i*h));
	  u[i][0] = sqrt(r)*sin(acos(-1 / r) / 2.0);
	 u[i][grid_size_x-1] = sqrt(r)*sin(acos(1 / r) / 2.0);
    }
}

void computeRes(std::vector<std::vector<double > > & u,
std::vector<std::vector<double > > & res_h,
std::vector<std::vector<double > > & f) {

  std::size_t grid_size_y = u.size();
  std::size_t grid_size_x = u[0].size();
   double h;

  h = 2.0/(double)(grid_size_x-1);
#pragma omp parallel for
  for(std::size_t i=1; i<grid_size_y; ++i){
    #pragma omp parallel for
	for(std::size_t j=1; j<grid_size_x-1; ++j){
	   if(!(i==(grid_size_y-1))){
		res_h[i][j] = f[i][j] - 1.0/pow(h,2)*(4.0*u[i][j] - u[i+1][j] -
		u[i-1][j] - u[i][j+1] - u[i][j-1]);
	   }
	   else if(j<(grid_size_x-1)/2){
	       res_h[i][j] = f[i][j] - 1.0/pow(h,2)*(4.0*u[i][j] - u[i-1][j] -
		u[i-1][j] - u[i][j+1] - u[i][j-1]);
	}

  }

}
}

void rbgs(std::vector<std::vector<double > > &u, int v1, std::vector<std::vector<double > > &f) {
//int num =omp_get_num_procs()
  std::size_t grid_size_y = u.size(); 
  std::size_t grid_size_x = u[0].size();

  double h = 2.0/(double)(grid_size_x-1);


  for(int it=0; it<v1; ++it){

    {

#pragma omp parallel for
	for(std::size_t i=1; i<grid_size_y; ++i){
	  #pragma omp parallel for
		for (std::size_t j=1;j<grid_size_x-1;++j){
			if((i+j)%2==0) {
			  if(!(i==grid_size_y-1)){
				double un=0.25*(u[i+1][j]+u[i-1][j]+u[i][j+1]+u[i][j-1]
				+ pow(h,2)*f[i][j]);
				u[i][j]=un;
			//	std::cout<<"thread number "<<omp_get_thread_num()<<std::endl;
			  }
			  else if(j<(grid_size_x-1)/2) {
			      double un=0.25*(u[i-1][j]+u[i-1][j]+u[i][j+1]+u[i][j-1]
				+ pow(h,2)*f[i][j]);
				u[i][j]=un;
			//	std::cout<<"thread number "<<omp_get_thread_num()<<std::endl;
			  }
			}
		  }
	  }
	  #pragma omp parallel for
	for(std::size_t i=1; i<grid_size_y; ++i){
	  #pragma omp parallel for
		for (std::size_t j=1;j<grid_size_x-1;++j){
			if((i+j)%2!=0) {
			  if(!(i==grid_size_y-1)){
				double un=0.25*(u[i+1][j]+u[i-1][j]+u[i][j+1]+u[i][j-1]
				+ pow(h,2)*f[i][j]);
				u[i][j]=un;
			//	std::cout<<"thread number "<<omp_get_thread_num()<<std::endl;
			  }
			  else if(j<(grid_size_x-1)/2) {
			      double un=0.25*(u[i-1][j]+u[i-1][j]+u[i][j+1]+u[i][j-1]
				+ pow(h,2)*f[i][j]);
				u[i][j]=un;
			//	std::cout<<"thread number "<<omp_get_thread_num()<<std::endl;
			  }
			}
		  }
	  }

    }
    
  }
}

void restriction(std::vector<std::vector<double > > & res_h,
std::vector<std::vector<double > >&f_2h){


  std::size_t coarser_grid_size_y = f_2h.size();
  std::size_t coarser_grid_size_x = f_2h[0].size();


  {
#pragma omp parallel for
  for(std::size_t i=0; i<coarser_grid_size_y; ++i){
    //  std::cout<<"thread number "<< omp_get_thread_num() <<std::endl;
    #pragma omp parallel for
    for(std::size_t j=0; j<coarser_grid_size_x; ++j){
      if(!(i==coarser_grid_size_y-1 && j>=(coarser_grid_size_x-1)/2)){

      if(i==0&&j==0){
	f_2h[i][j] = 1.0/9.0*(res_h[2*i+1][2*j+1]+2*res_h[2*i+1][2*j]+4*res_h[2*i][2*j]+2*res_h[2*i][2*j+1]);
      }
      else if(i==0&&j==(coarser_grid_size_x-1)){
	f_2h[i][j] = 1.0/9.0*(2*res_h[2*i+1][2*j]+res_h[2*i+1][2*j-1]+2*res_h[2*i][2*j-1]+4*res_h[2*i][2*j]);
      }
      else if(i==(coarser_grid_size_y-1)&&j==0){
	f_2h[i][j] = 1.0/12.0*(res_h[2*i-1][2*j+1]+2*res_h[2*i-1][2*j]+4*res_h[2*i][2*j]+2*res_h[2*i][2*j+1]+2*res_h[2*i-1][2*j]+res_h[2*i-1][2*j+1]);
      }
      else if(i==(coarser_grid_size_y-1)&&j==(coarser_grid_size_x-1)){
	f_2h[i][j] = 1.0/12.0*(2*res_h[2*i-1][2*j]+res_h[2*i-1][2*j-1]+2*res_h[2*i][2*j-1]+4*res_h[2*i][2*j]+res_h[2*i-1][2*j-1]+2*res_h[2*i-1][2*j]);
      }
      else if(i==0 && j!=0 && j!=(coarser_grid_size_x-1)){
	f_2h[i][j] = 1.0/12.0*(res_h[2*i+1][2*j+1]+2*res_h[2*i+1][2*j]+res_h[2*i+1][2*j-1]+2*res_h[2*i][2*j-1]+4*res_h[2*i][2*j]+2*res_h[2*i][2*j+1]);
      }
      else if(i==(coarser_grid_size_y-1) && j!=0 && j!=(coarser_grid_size_x-1)){
	f_2h[i][j] = 1.0/16.0*(res_h[2*i-1][2*j+1]+2*res_h[2*i-1][2*j]+res_h[2*i-1][2*j-1]+2*res_h[2*i][2*j-1]+4*res_h[2*i][2*j]+2*res_h[2*i][2*j+1]+res_h[2*i-1][2*j-1]+2*res_h[2*i-1][2*j]+res_h[2*i-1][2*j+1]);
      }
      else if(j==0 && i!=0 && i!=(coarser_grid_size_y-1)){
	f_2h[i][j] = 1.0/12.0*(res_h[2*i+1][2*j+1]+2*res_h[2*i+1][2*j]+4*res_h[2*i][2*j]+2*res_h[2*i][2*j+1]+2*res_h[2*i-1][2*j]+res_h[2*i-1][2*j+1]);
      }
      else if(j==(coarser_grid_size_x-1) && i!=0 && i!=(coarser_grid_size_y-1)){
	f_2h[i][j] = 1.0/12.0*(2*res_h[2*i+1][2*j]+res_h[2*i+1][2*j-1]+2*res_h[2*i][2*j-1]+4*res_h[2*i][2*j]+res_h[2*i-1][2*j-1]+2*res_h[2*i-1][2*j]);
      }
      else {
	f_2h[i][j] = 1.0/16.0*(res_h[2*i+1][2*j+1]+2*res_h[2*i+1][2*j]+res_h[2*i+1][2*j-1]+2*res_h[2*i][2*j-1]+4*res_h[2*i][2*j]+2*res_h[2*i][2*j+1]+res_h[2*i-1][2*j-1]+2*res_h[2*i-1][2*j]+res_h[2*i-1][2*j+1]);
      }
      }

     }
  }
  }
    
  }


void correction(std::vector<std::vector<double > >& u,
std::vector<std::vector<double > >&c_2h){


  std::size_t coarser_grid_size_y = c_2h.size();
  std::size_t coarser_grid_size_x = c_2h[0].size();
  std::size_t grid_size_y = u.size();
  std::size_t grid_size_x = u[0].size();
#pragma omp parallel for
  for (std::size_t i=0; i<coarser_grid_size_y; ++i) {
    #pragma omp parallel for
    for (std::size_t j = 0;j<coarser_grid_size_x-1; ++j) {
	u[2*i][2*j] += c_2h[i][j];
      	if((2*i+1)<grid_size_y && 2*j!=0 && 2*j!=grid_size_x-1) {
	  u[2*i+1][2*j] += 1.0/2.0*(c_2h[i][j] + c_2h[i+1][j]);
      	}
      	if((2*j+1)<grid_size_x-1&& 2*i!=0 && 2*i!= grid_size_x-1){
	  u[2*i][2*j+1] += 1.0/2.0*(c_2h[i][j] + c_2h[i][j+1]);
      	}
      	if((2*i+1)<grid_size_y && (2*j+1)<grid_size_x-1) {
	  u[2*i+1][2*j+1] += 1.0/4.0*(c_2h[i][j]+c_2h[i+1][j] +
			      c_2h[i][j+1] +c_2h[i+1][j+1]);
      	}
   }
  }
}


void interpolation (std::vector<std::vector<double > >& u,
std::vector<std::vector<double > >&c_2h){

  std::size_t coarser_grid_size_y = c_2h.size();
  std::size_t coarser_grid_size_x = c_2h[0].size();
  std::size_t grid_size_y = u.size();
  std::size_t grid_size_x = u[0].size();
#pragma omp parallel for
  for (std::size_t i=0; i<coarser_grid_size_y; ++i) {
    #pragma omp parallel for
    for (std::size_t j = 0;j<coarser_grid_size_x-1; ++j) {
	u[2*i][2*j] = c_2h[i][j];
      	if((2*i+1)<grid_size_y && 2*j!=0 && 2*j!=grid_size_x-1) {
	  u[2*i+1][2*j] = 1.0/2.0*(c_2h[i][j] + c_2h[i+1][j]);
      	}
      	if((2*j+1)<grid_size_x-1&& 2*i!=0 && 2*i!= grid_size_x-1){
	  u[2*i][2*j+1] = 1.0/2.0*(c_2h[i][j] + c_2h[i][j+1]);
      	}
      	if((2*i+1)<grid_size_y && (2*j+1)<grid_size_x-1) {
	  u[2*i+1][2*j+1] = 1.0/4.0*(c_2h[i][j]+c_2h[i+1][j] +
			      c_2h[i][j+1] +c_2h[i+1][j+1]);
      	}
   }
  }
}


