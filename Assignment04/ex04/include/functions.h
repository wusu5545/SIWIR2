#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <boost/make_shared.hpp>
#include <boost/shared_ptr.hpp>

using namespace std;

void Initialization(string parfile,string datfile);

template<typename Type>
Type operator[](size_t x,size_t y) const
{
  return (*this)[3*x+y];
}

template<typename Type>
Type & operator[](size_t x,size_t y)
{
  return (*this)[3*x+y];
}


#endif//FUNCTIONS_H
