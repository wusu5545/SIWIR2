#ifndef FUNCTIONS
#define FUNCTIONS

#include "GRID.h"
using namespace lbm;

void vtk_output(int vtk_out_freq, int iteration, std::string & vtk_out_name_base, int& vtk_out_count,
std::size_t xsize,std::size_t ysize, const FLAGS& flags,const D_Field& d, const V_Field& v);

#endif