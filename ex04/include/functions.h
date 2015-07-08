#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <boost/make_shared.hpp>
#include <boost/shared_ptr.hpp>

#define x(i,j) x[i*3+j]
#define v(i,j) v[i*3+j]
#define F(i,j) F[i*3+j]
#define F_old(i,j) F_old[i*3+j]
#define r(i,j) r[i*3+j]
#define cell_num(i,j) cell_num[i*3+j]

using namespace std;

void LinkedCell();
void ForceCalculation();
void Position_update();
void Velocity_update();
void vtk_output(size_t iteration,size_t &vtk_count);
void Initialization(string parfile,string datfile);

#endif//FUNCTIONS_H
