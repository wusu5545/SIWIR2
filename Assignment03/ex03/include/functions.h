#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include <fstream>
#include <iostream>
#include "filereader.h"
#include "lbm.h"

using namespace std;
using namespace lbm;

void vtk_output(size_t vtk_step,size_t t,string &vtk_file,size_t &vtk_count,
		size_t sizex,size_t sizey,const Flags &flags,const D_Field& rho,const V_Field &v);
void Solver(FileReader &reader);

#endif//FUNCTIONS_H
