#include <assert.h>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <math.h>
#include <omp.h>

#include "stencil.h"
#include "mgsolver.h"
#include "grid.h"
#include "util.h"
#include "config.h"
#include "cgsolver.h"


#ifdef MEASURE_TIME_PER_FUNCTION
#include "timer.h"
#endif


MGSolver::MGSolver() {
	real factor = 1/real(16);

	Stencil restrictionStencil{
		4 * factor,
		  2 * factor,
		  2 * factor,
		  2 * factor,
		  2 * factor,
		  1 * factor,
		  1 * factor,
		  1 * factor,
		  1 * factor
	};

	setRestrictionStencil(restrictionStencil);

	factor = 1/real(4);
	Stencil prolongationStencil{
		4 * factor,
		  2 * factor,
		  2 * factor,
		  2 * factor,
		  2 * factor,
		  1 * factor,
		  1 * factor,
		  1 * factor,
		  1 * factor
	};

	setProlongationStencil(prolongationStencil);
}

MGSolver::~MGSolver() {
	//nothing to do
}

#ifdef MEASURE_TIME_PER_FUNCTION

double t_smooth = 0;
double t_calcResidual = 0; 
double t_restrictGrid = 0;
double t_prolongateAndCorrect = 0;
double t_cgsolve = 0;

#endif

void MGSolver::solve (Grids& x_, Grids& f_, Grids& r_, Grid& d, Grid& z, size_t levels, size_t pre, size_t post) {
	//CGSolver::solve(x_[levels-1], f_[levels-1], r_[levels-1], d, z, CGSTEPS);
	//if(d.hsize() != -1) return;
	
	if(d.hsize() == -1){

		//Smoother::rb_gauss_seidel_RBAtOnce(x, f, 1, 0);
		Smoother::rb_gauss_seidel_f(x_[levels-1], f_[levels-1], 1, 0);

		return;
	}



#if MEASUREL2NORM == 1
	Grid r(levels);
	calcResidual_f(x_[levels-1],f_[levels-1],r);
	double newL2Norm, oldL2Norm = r.calcL2Norm();
#endif



	/************ VCYCLES *****************/		
	for(size_t cycle = 0; cycle < FLOATCYCLES ; ++cycle){
		// --- downward phase

		for(size_t level = levels-1; level > COARSEST; level--) {
			smooth_f(x_[level], f_[level], pre, cycle);
			calcResidual_f(x_[level], f_[level], r_[level]);
			restrictGrid_f(r_[level], f_[level-1]);
			x_[level-1].fill_f(0);
		}

		// --- coarsest grid ---
#if MEASURE_TIME_PER_FUNCTION == 1
		siwir::Timer timer;
#endif
#if COARSEST_SOLVER == 0
		smooth_f(x_[COARSEST], f_[COARSEST], 1, cycle);
#elif COARSEST_SOLVER == 1
		CGSolver::solve(x_[COARSEST], f_[COARSEST], r_[COARSEST], d, z, CGSTEPS);
#endif
#if MEASURE_TIME_PER_FUNCTION == 1
		t_cgsolve += timer.elapsed();
#endif

		// --- upward phase ---
		for(size_t level = COARSEST+1; level < levels; ++level){			
			prolongateAndCorrect_f(x_[level-1], x_[level]);
			smooth_f(x_[level], f_[level], post, cycle);
		}


#if MEASUREL2NORM == 1
		calcResidual_f(x_[levels-1], f_[levels-1], r_[levels-1]);
		newL2Norm = r_[levels-1].calcL2Norm_f(); 
		std::cout << " Cycle " << std::setw(3) << cycle+1 << ": "  
			<< "L2-norm of the residual = " << std::setw(11) << newL2Norm << " "
			<< "Convergence rate = " << newL2Norm / oldL2Norm << std::endl;
		oldL2Norm = newL2Norm;

		if(newL2Norm < 0.0000918){
			break;
		}
#endif
	}



#if MEASURE_TIME_PER_FUNCTION == 1
	std::cout << std::setw(20) << "smooth: " << std::setw(11) << t_smooth << std::endl;
	std::cout << std::setw(20) << "residual: " << std::setw(11) << t_calcResidual << std::endl;
	std::cout << std::setw(20) << "restrict: " << std::setw(11) << t_restrictGrid << std::endl;
	std::cout << std::setw(20) << "prolongate: " << std::setw(11) << t_prolongateAndCorrect << std::endl;
	std::cout << std::setw(20) << "solve coarsest: " << std::setw(11) << t_cgsolve << std::endl;
#endif
}

void MGSolver::setSmootherStencil (const Stencil& stencil) {
	//simple struct copy is handled automatically in c++
	A_ = stencil;
}
void MGSolver::setRestrictionStencil (const Stencil& stencil) {
	//simple struct copy is handled automatically in c++
	I_ = stencil;
}
void MGSolver::setProlongationStencil (const Stencil& stencil) {
	//simple struct copy is handled automatically in c++
	J_ = stencil;
}


void MGSolver::smooth_f (Grid& x, Grid& f, size_t iter, size_t cycle) {
#if MEASURE_TIME_PER_FUNCTION == 1
	siwir::Timer timer;
#endif

	Smoother::rb_gauss_seidel_f(x, f, iter, cycle);

#if MEASURE_TIME_PER_FUNCTION == 1
	t_smooth += timer.elapsed();
#endif
}

void MGSolver::calcResidual_f(const Grid& x, const Grid& f, Grid& r) {
#if MEASURE_TIME_PER_FUNCTION == 1
	siwir::Timer timer;
#endif

	const int halfSize = x.getHalfSize();

	const float H = x.hsize();
	const float center = real(4) / (H * H);
	const float outer  = real(1) / (H * H);

#pragma omp parallel
	// i < 0
#pragma omp for nowait
	for (int i = -halfSize+1; i < 0; ++i) {
		for (int j = -halfSize+1; j < halfSize; ++j) {
			r(i, j) = f(i, j) - center * x(i,j) 
				+ outer * (x(i,j+1) + x(i,j-1) + 
						x(i+1,j) + x(i-1,j)); 
		}
	}
	// i = 0
#pragma omp for nowait
	for (int j = -halfSize+1; j < 0; ++j) {
		r(0, j) = f(0, j) - center * x(0,j) 
			+ outer * (x(0,j+1) + x(0,j-1) + 
					x(-1,j) + x(-1,j)); 
	}

#if MEASURE_TIME_PER_FUNCTION == 1
	t_calcResidual += timer.elapsed();
#endif
}

void MGSolver::restrictGrid_f(const Grid& r, Grid& f) {
#if MEASURE_TIME_PER_FUNCTION == 1
	siwir::Timer timer;
#endif

	const int halfSize = f.getHalfSize();

#pragma omp parallel for
	for(int i = -halfSize+1; i < 0; ++i){
		for(int j = -halfSize+1; j < halfSize; ++j){
			f(i,j) =
				I_.C  *	r(i*2	, j*2	) +
				I_.N  *	r(i*2+1	, j*2	) +
				I_.S  *	r(i*2-1	, j*2	) +
				I_.W  *	r(i*2	, j*2-1	) +
				I_.E  *	r(i*2	, j*2+1	) +
				I_.NW * r(i*2+1	, j*2-1	) +
				I_.NE *	r(i*2+1	, j*2+1	) +
				I_.SW *	r(i*2-1	, j*2-1	) +
				I_.SE *	r(i*2-1	, j*2+1	) ; 
		}
	}

	for(int j = -halfSize+1; j < 0; ++j){
		f(0,j) =
			I_.C  *	r(0 , j*2	) +
			I_.N  *	r(-1, j*2	) +
			I_.S  *	r(-1, j*2	) +
			I_.W  *	r(0 , j*2-1	) +
			I_.E  *	r(0 , j*2+1	) +
			I_.NW * r(-1, j*2-1	) +
			I_.NE *	r(-1, j*2+1	) +
			I_.SW *	r(-1, j*2-1	) +
			I_.SE *	r(-1, j*2+1	) ; 
	}

#if MEASURE_TIME_PER_FUNCTION == 1
	t_restrictGrid += timer.elapsed();
#endif
}


void MGSolver::calcResidual_d(const Grid& x, const Grid& f, Grid& r) {
#if MEASURE_TIME_PER_FUNCTION == 1
	siwir::Timer timer;
#endif

	const int halfSize = x.getHalfSize();

	const double H = x.hsize();
	const double center = 4.0 / (H * H);
	const double outer  = 1.0 / (H * H);

#pragma omp parallel
{
	if(halfSize > MINHALFSIZE){
		// i < 0
		int jj = -halfSize+1;
		for(; jj < halfSize-BLOCKSIZE_RESIDUAL; jj+=BLOCKSIZE_RESIDUAL){
#pragma omp for nowait
			for (int i = -halfSize+1; i < 0; ++i) {
				for (int j = jj; j < jj+BLOCKSIZE_RESIDUAL; ++j) {
					r(i, j, 'd') = f(i, j, 'd') - center * x(i,j, 'd') 
						+ outer * (x(i,j+1, 'd') + x(i,j-1,'d') + 
								x(i+1,j,'d') + x(i-1,j,'d'));;
				}
			}
		}
	
#pragma omp for nowait
		for (int i = -halfSize+1; i < 0; ++i) {
			for (int j = jj; j < halfSize; ++j) {
				r(i, j, 'd') = f(i, j, 'd') - center * x(i,j, 'd') 
					+ outer * (x(i,j+1, 'd') + x(i,j-1,'d') + 
							x(i+1,j,'d') + x(i-1,j,'d'));;
			}
		}
	} else {
#pragma omp for nowait
		for (int i = -halfSize+1; i < 0; ++i) {
			for (int j = -halfSize+1; j < halfSize; ++j) {
				r(i, j, 'd') = f(i, j, 'd') - center * x(i,j, 'd') 
					+ outer * (x(i,j+1, 'd') + x(i,j-1,'d') + 
							x(i+1,j,'d') + x(i-1,j,'d'));;
			}
		}
	}

	// i = 0
#pragma omp for nowait
	for (int j = -halfSize+1; j < 0; ++j) {
		r(0, j,'d') = f(0, j,'d') - center * x(0,j,'d') 
			+ outer * (x(0,j+1,'d') + x(0,j-1,'d') + 
					x(-1,j,'d') + x(-1,j,'d')); 
	}
}
#if MEASURE_TIME_PER_FUNCTION == 1
	t_calcResidual += timer.elapsed();
#endif
}


void MGSolver::calcResidual_fTod(const Grid& x, const Grid& f, Grid& r) {
#if MEASURE_TIME_PER_FUNCTION == 1
	siwir::Timer timer;
#endif

	const int halfSize = x.getHalfSize();

	const double H = x.hsize();
	const double center = 4.0 / (H * H);
	const double outer  = 1.0 / (H * H);

#pragma omp parallel
	// i < 0
#pragma omp for nowait
	for (int i = -halfSize+1; i < 0; ++i) {
		for (int j = -halfSize+1; j < halfSize; ++j) {
			r(i, j, 'd') = f(i, j) - center * x(i,j, 'd') 
				+ outer * (x(i,j+1, 'd') + x(i,j-1,'d') + 
						x(i+1,j,'d') + x(i-1,j,'d')); 
		}
	}
	// i = 0
#pragma omp for nowait
	for (int j = -halfSize+1; j < 0; ++j) {
		r(0, j,'d') = f(0, j) - center * x(0,j,'d') 
			+ outer * (x(0,j+1,'d') + x(0,j-1,'d') + 
					x(-1,j,'d') + x(-1,j,'d')); 
	}

#if MEASURE_TIME_PER_FUNCTION == 1
	t_calcResidual += timer.elapsed();
#endif
}



void MGSolver::restrictGrid_d(const Grid& r, Grid& f) {
#if MEASURE_TIME_PER_FUNCTION == 1
	siwir::Timer timer;
#endif

	const int halfSize = f.getHalfSize();

#pragma omp parallel for
	for(int i = -halfSize+1; i < 0; ++i){
		for(int j = -halfSize+1; j < halfSize; ++j){
			f(i,j,'d') =
				I_.C  *	r(i*2	, j*2	,'d') +
				I_.N  *	r(i*2+1	, j*2	,'d') +
				I_.S  *	r(i*2-1	, j*2	,'d') +
				I_.W  *	r(i*2	, j*2-1	,'d') +
				I_.E  *	r(i*2	, j*2+1	,'d') +
				I_.NW * r(i*2+1	, j*2-1	,'d') +
				I_.NE *	r(i*2+1	, j*2+1	,'d') +
				I_.SW *	r(i*2-1	, j*2-1	,'d') +
				I_.SE *	r(i*2-1	, j*2+1	,'d') ; 
		}
	}

	for(int j = -halfSize+1; j < 0; ++j){
		f(0,j,'d') =
			I_.C  *	r(0 , j*2	,'d') +
			I_.N  *	r(-1, j*2	,'d') +
			I_.S  *	r(-1, j*2	,'d') +
			I_.W  *	r(0 , j*2-1	,'d') +
			I_.E  *	r(0 , j*2+1	,'d') +
			I_.NW * r(-1, j*2-1	,'d') +
			I_.NE *	r(-1, j*2+1	,'d') +
			I_.SW *	r(-1, j*2-1	,'d') +
			I_.SE *	r(-1, j*2+1	,'d') ; 
	}

#if MEASURE_TIME_PER_FUNCTION == 1
	t_restrictGrid += timer.elapsed();
#endif
}


void MGSolver::prolongateAndCorrect_f(const Grid& e, Grid& x) {
#if MEASURE_TIME_PER_FUNCTION == 1
	siwir::Timer timer;
#endif

	const int halfSize = e.getHalfSize();


#pragma omp parallel for		
	for(int i = -halfSize+1; i < 0; ++i){
		for(int j = -halfSize+1; j < halfSize; ++j){
			x(i<<1    , j<<1	) 	+= J_.C  *  e(i	 , j);

			x((i<<1)  , (j<<1)+1)	+= J_.E  * (e(i	 , j  ) +
					e(i  , j+1));

			x((i<<1)+1, (j<<1)  )	+= J_.N  * (e(i  , j  ) +
					e(i+1, j  ));

			x((i<<1)+1, (j<<1)+1)   += J_.NE * (e(i  , j  ) +
					e(i+1, j  ) +
					e(i  , j+1) +
					e(i+1, j+1));
		}	
	}

	for(int j = -halfSize+1; j < 0; ++j){
		x(0, j<<1	) 	+= J_.C  * e(0, j);

		x(0, (j<<1)+1)	+= J_.E  * (e(0, j) + e(0, j+1));
	}

#if MEASURE_TIME_PER_FUNCTION == 1
	t_prolongateAndCorrect += timer.elapsed();
#endif
}

void MGSolver::prolongateAndCorrect_d(const Grid& e, Grid& x) {
#if MEASURE_TIME_PER_FUNCTION == 1
	siwir::Timer timer;
#endif

	const int halfSize = e.getHalfSize();


#pragma omp parallel for		
	for(int i = -halfSize+1; i < 0; ++i){
		for(int j = -halfSize+1; j < halfSize; ++j){
			x(i<<1    , j<<1	,'d') 	+= J_.C  *  e(i	 , j,'d');

			x((i<<1)  , (j<<1)+1,'d')	+= J_.E  * (e(i	 , j  ,'d') + e(i  , j+1,'d'));

			x((i<<1)+1, (j<<1)  ,'d')	+= J_.N  * (e(i  , j  ,'d') + e(i+1, j  ,'d'));

			x((i<<1)+1, (j<<1)+1,'d')   += J_.NE * (e(i  , j  ,'d') +
					e(i+1, j  ,'d') +
					e(i  , j+1,'d') +
					e(i+1, j+1,'d'));
		}	
	}

	for(int j = -halfSize+1; j < 0; ++j){
		x(0, j<<1	,'d') 	+= J_.C  * e(0, j,'d');

		x(0, (j<<1)+1,'d')	+= J_.E  * (e(0, j,'d') + e(0, j+1,'d'));
	}

#if MEASURE_TIME_PER_FUNCTION == 1
	t_prolongateAndCorrect += timer.elapsed();
#endif
}

