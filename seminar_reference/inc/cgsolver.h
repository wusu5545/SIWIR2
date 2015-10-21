#ifndef CGSOLVER_H
#define CGSOLVER_H

#include "grid.h"

class CGSolver {
	public:
	static void solve(Grid& x, Grid& f, Grid& r, Grid& d, Grid& z, size_t iter);
};

#endif
