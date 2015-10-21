#ifndef SMOOTHER_H 
#define SMOOTHER_H

#include "grid.h"

class Smoother {
	public:
		static void rb_gauss_seidel_f(Grid& a, const Grid& f, size_t iter, size_t cycle);
		static void rb_gauss_seidel_f_RBAtOnce(Grid& a, const Grid& f, size_t iter, size_t cycle);
};

#endif
