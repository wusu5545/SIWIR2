#include "cgsolver.h"
#include "grid.h"
#include "mgsolver.h"
#include "types.h"


void CGSolver::solve(Grid& x, Grid& f, Grid& r, Grid& d, Grid& z, size_t iter){
#if USE_CG == 1	
	const int halfSize = r.getHalfSize();
	const real H = x.hsize();
	const real center = real(4) / (H * H);
	const real outer  = real(1) / (H * H);

	MGSolver::calcResidual(x, f, r);
	for(int i = -halfSize; i <= 0; ++i){
		for(int j = -halfSize; j < halfSize; ++j){
			d(i, j) = r(i, j);
		}
	}
	
	real tmp = -1;

	for(size_t k = 0; k < iter; ++k){
		
		// --- z = A*d ---	
		// i < 0
#pragma omp parallel for
		for (int i = -halfSize+1; i < 0; ++i) {
			for (int j = -halfSize+1; j < halfSize; ++j) {
				z(i, j) = center * d(i,j) 
						- outer * (d(i,j+1) + d(i,j-1) + d(i+1,j) + d(i-1,j)); 
			}
		}
		// i = 0
#pragma omp parallel for
		for (int j = -halfSize+1; j < 0; ++j) {
			z(0, j) = center * d(0,j) 
					- outer * (d(0,j+1) + d(0,j-1) + d(-1,j) + d(-1,j)); 
		}
		
		// --- alphak berechnen ---

		//if(tmp == -1){
			tmp = 0;
			for (int i = -halfSize; i < 0; ++i) {
				for (int j = -halfSize; j < halfSize+1; ++j) {
					tmp += r(i, j) * r(i, j);
				}
			}
			
			for (int j = -halfSize; j < 0; ++j) {
				tmp += r(0, j) * r(0, j);
			}
			
			for (int i = 1; i < halfSize+1; ++i) {
				for (int j = -halfSize; j < halfSize+1; ++j) {
					tmp += r(-i, j) * r(-i, j);
				}
			}
		//}
		real tmp2 = 0;
		for (int i = -halfSize; i < 0; ++i) {
			for (int j = -halfSize; j < halfSize+1; ++j) {
				tmp2 += z(i, j) * d(i, j);
			}
		}
		
		for (int j = -halfSize; j < 0; ++j) {
			tmp2 += z(0, j) * d(0, j);
		}
		
		for (int i = 1; i < halfSize+1; ++i) {
			for (int j = -halfSize; j < halfSize+1; ++j) {
				tmp2 += z(-i, j) * d(-i, j);
			}
		}
		real alphak = tmp/tmp2;

		//std::cout << alphak << " " << tmp << " " << tmp2 << std::endl;

		// --- x = x + alphak*d
		for (int i = -halfSize; i < 0; ++i) {
			for (int j = -halfSize; j < halfSize+1; ++j) {
				x(i, j) += alphak*d(i, j);
			}
		}
		
		for (int j = -halfSize; j < 0; ++j) {
			x(0, j) += alphak*d(0, j);
		}

		// --- r = r - alphak*z
		for (int i = -halfSize; i < 0; ++i) {
			for (int j = -halfSize; j < halfSize+1; ++j) {
				r(i, j) -= alphak*z(i, j);
			}
		}
		
		for (int j = -halfSize; j < 0; ++j) {
			r(0, j) -= alphak*z(0, j);
		}

		// --- betak berechnen ---
		tmp2 = 0;
		for (int i = -halfSize; i < 0; ++i) {
			for (int j = -halfSize; j < halfSize+1; ++j) {
				tmp2 += r(i, j) * r(i, j);
			}
		}
	
		for (int j = -halfSize; j < 0; ++j) {
			tmp2 += r(0, j) * r(0, j);
		}
		
		for (int i = 1; i < halfSize+1; ++i) {
			for (int j = -halfSize; j < halfSize+1; ++j) {
				tmp2 += r(-i, j) * r(-i, j);
			}
		}
		real betak = tmp2/tmp;
		tmp = tmp2;

		// --- d = r + betak*d ---
		for (int i = -halfSize; i < 0; ++i) {
			for (int j = -halfSize; j < halfSize+1; ++j) {
				d(i, j) = r(i, j) + betak*d(i, j);
			}
		}
		for (int j = -halfSize; j < 0; ++j) {
			d(0, j) = r(0, j) + betak*d(0, j);
		}
	}
#endif
}
