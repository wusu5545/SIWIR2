#ifndef FUNCTIONS
#define FUNCTIONS

#include <vector>
#include <fstream>

void rbgs(std::vector<std::vector<double > >& u, int v1,
std::vector<std::vector<double > > &f);
void computeRes(std::vector<std::vector<double > > & u,
std::vector<std::vector<double > > & res_h,
std::vector<std::vector<double > > & f);
void restriction(std::vector<std::vector<double > >& res_h,
std::vector<std::vector<double > >& f_2h);
void mgm(std::vector<std::vector<double > > &u, std::vector<std::vector<double > > &f, int v1, int v2, int gamma,
	 int l, std::size_t grid_size);
void correction(std::vector<std::vector<double > >& u,
std::vector<std::vector<double > >&c_2h);
void init_and_boundary(std::vector<std::vector <double> > & u);
void solveMG(int l, std::vector<std::vector<double > >& u, std::vector<std::vector<double > >& f);
void fmg(std::vector<std::vector<double > > &u, std::vector<std::vector<double > > &f,int v1, int v2, int gamma, int miu,
	 int l, std::size_t grid_size);
void interpolation (std::vector<std::vector<double > >& u,
std::vector<std::vector<double > >&c_2h);

#endif
