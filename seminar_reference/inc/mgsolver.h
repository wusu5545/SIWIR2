#ifndef MGSOLVER_H
#define MGSOLVER_H

#include <vector>

#include "smoother.h"
#include "stencil.h"
#include "grid.h"



class MGSolver {
	private:
		//Type definitions
		
	public:
		
		//Constructors
		explicit MGSolver();
		
		//Destructor
		~MGSolver();
		
		//Multigrid functions
		void solve (Grids& x_, Grids& f_, Grids& r_, Grid& d, Grid& z, size_t l, size_t pre, size_t post);
		
		//Utility functions
		void setSmootherStencil (const Stencil& stencil);
		void setRestrictionStencil (const Stencil& stencil);
		void setProlongationStencil (const Stencil& stencil);
		
		static void calcResidual_f(const Grid& x, const Grid& f, Grid& r);
		static void calcResidual_d(const Grid& x, const Grid& f, Grid& r);
		static void calcResidual_fTod(const Grid& x, const Grid& f, Grid& r);
	private:
		//Multigrid functions
		void smooth_f (Grid& x, Grid& f, size_t iter, size_t cycle);
		void restrictGrid_f(const Grid& r, Grid& f);
		void restrictGrid_d(const Grid& r, Grid& f);
		void prolongateAndCorrect_f(const Grid& e, Grid& x);
		void prolongateAndCorrect_d(const Grid& e, Grid& x);
		
		//Member variables
		
		Stencil A_; //The smoothing stencil.
		Stencil I_; //The restriction stencil.
		Stencil J_; //The prolongation stencil.
};

#endif
