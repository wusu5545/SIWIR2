
#ifndef PRECISION_H
#define PRECISION_H
#include <config.h>
//defines the precision of real numbers


#if DOUBLE_PRECISION == 1
typedef double real;
#else
typedef float real;
#endif

#endif
