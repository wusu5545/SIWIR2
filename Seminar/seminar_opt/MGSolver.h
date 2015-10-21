#ifndef MGSOLVER_H
#define MGSOLVER_H

#include<boost/function.hpp>
#include<iostream>
#include<iomanip> 
#include<assert.h>
#include<boost/make_shared.hpp>
#include<boost/shared_ptr.hpp>
#include<math.h>
#include<omp.h>

using namespace std;
using boost::make_shared;
using boost::shared_ptr;
using namespace boost;

/* !\ brief Floating point data type for the multigrid assignment .
 //
 // This type definition offers the possibility to switch the floating
 // point precision of the multigrid assignment between float , double
 // and long double .
 //
 // Valid types for the \a real floating point type :
 // <a>float , double , long double </a>
 */
typedef double real;

struct Stencil
{
   // Constructor
   explicit inline Stencil();
   inline Stencil(real c,real n,real s,real w,real e,real nw,real ne,real sw,real se);
   // Member variables
   real C; // The central stencil entry .
   real N; // The stencil entry for access to the northern neighbor .
   real S; // The stencil entry for access to the southern neighbor .
   real W; // The stencil entry for access to the western neighbor .
   real E; // The stencil entry for access to the eastern neighbor .
   real NW; // The stencil entry for access to the north - western neighbor .
   real NE; // The stencil entry for access to the north - eastern neighbor .
   real SW; // The stencil entry for access to the south - western neighbor .
   real SE; // The stencil entry for access to the south - eastern neighbor .
};
inline Stencil::Stencil()
{
  Stencil(1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1);
}
inline Stencil::Stencil(real c, real n, real s, real w, real e, real nw, real ne, real sw, real se)
{
    this->C = c;
    this->N = n;
    this->S = s;
    this->W = w;
    this->E = e;
    this->NW = nw;
    this->NE = ne;
    this->SW = sw;
    this->SE = se;
}

class Grid
{

public:
   // Constructors
   explicit Grid();
   explicit Grid( size_t level );
   Grid( const Grid & grid );

   // Destructor
   ~ Grid();

   // Copy assignment
// 	const Grid & operator =(const Grid & grid);

   // Get functions
   inline size_t xsize() const;
   inline size_t ysize() const;
   inline real hsize() const;

   // Arithmetic operators
// 	Grid & operator +=(const Grid & rhs);
// 	Grid & operator -=(const Grid & rhs);
   // Access functions
   inline real & operator()( size_t i, size_t j );
   inline real operator()( size_t i, size_t j ) const;

   // Utility functions
   real calcL2Norm() const;
   void fill( real value );
   void setBoundary( real value );
   void resize( size_t level );
   void resize( size_t y, size_t x );

private:
   // Member variables
   bool CASE;
   size_t y_; // The current number of rows of the 2D grid
   size_t x_; // The current number of columns of the 2D grid
   real h_; // The mesh size of the grid
   boost::shared_ptr< real[] > v_; // The values contained in the grid
};

Grid::Grid()
{
   resize( 0 );
}

Grid::Grid( size_t level )
{
   resize( level );
}

Grid::Grid( const Grid& grid )
{
   ( *this ) = grid;
}

Grid::~Grid()
{

}

inline size_t Grid::xsize() const
{
   return x_;
}

inline size_t Grid::ysize() const
{
   return y_;
}

inline real Grid::hsize() const
{
   return h_;
}

real Grid::calcL2Norm() const
{
   real res = 0;
   for( size_t i = 0; i < y_; i++ )
      for( size_t j = 0; j < x_; j++ )
      {
	res += v_[i * x_ + j] * v_[i * x_ + j];
      }
   res = sqrt( res / ( x_ * y_ ) );

   return res;
}

void Grid::fill( real value )
{
   v_ = make_shared< real[] >( x_ * y_, value );
}

void Grid::setBoundary( real value )
{

}

void Grid::resize( size_t level )
{
   x_ = pow( 2, level ) + 1;
   y_ = pow( 2, level -1) + 1;
   v_ = make_shared< real[] >( x_ * y_ );
   h_ = real( 2.0 ) / ( pow( real( 2.0 ), level ) ); //h
}

void Grid::resize( size_t y, size_t x )
{
   x_ = x;
   y_ = y;
   v_ = make_shared< real[] >( x_ * y_ );
}

/* !\ brief Returns the current grid value at the specified position .
 //
 // \ param i The index in y- direction .
 // \ param j The index in x- direction .
 // \ return The grid value at the specified position .
 //
 // \b Note : In debug mode the two indices are checked via assertion !
 */
inline real& Grid::operator()( size_t i, size_t j )
{
//    assert( i < y_ );
//    assert( j < x_ );
   return v_[i * x_ + j];
}

inline real Grid::operator()( size_t i, size_t j ) const
{
//    assert( i < y_ );
//    assert( j < x_ );
   return v_[i * x_ + j];
}

class MGSolver
{
private:
   // Type definitions
   typedef boost::shared_ptr< Grid[] > Grids;

public:
   // Constructors
   explicit MGSolver( size_t nlevels );
   // Destructor
   ~ MGSolver();
   //Set number of vcycles
   void setNumberOfCycle();
   //Set Ghost layer
   void setGhostLayer( size_t addy, size_t addx );
   // Multigrid functions
   void solve( Grid & x, const Grid & f, size_t cycles, size_t pre, size_t post ,size_t levels);
   void solveFMG( Grid& x, const Grid& f, size_t pre, size_t post );
   // Utility functions
   void setSmootherStencil( const Stencil & stencil );
   void setRestrictionStencil( const Stencil & stencil );
   void setProlongationStencil( const Stencil & stencil );
private:
   boost::shared_ptr<size_t[]> N_;
   // Multigrid functions
   void smooth( Grid & x, const Grid & f, size_t iter );
   void calcResidual( const Grid & x, const Grid & f, Grid & r );
   void restrict( const Grid & r, Grid & f );
   void prolongateAndCorrect( const Grid & e, Grid & x );

   // Member variables
   size_t levels_; // The maximum number of levels .
   size_t vcycles_;	// The number of V-cycles
   Grids x_; // The grid hierarchy for the grids of unknowns .
   Grids f_; // The grid hierarchy for the right - hand side grids .
   Grids r_; // The grid hierarchy for the residual grids .
   Stencil A_; // The red - black Gauss - Seidel smoothing stencil .
   Stencil I_; // The restriction stencil .
   Stencil J_; // The prolongation stencil .
};

MGSolver::MGSolver( size_t nlevels )
{

   levels_ = nlevels;
   N_ = make_shared<size_t[]>(levels_);
   setNumberOfCycle();
   x_ = make_shared< Grid[] >( levels_ );
   f_ = make_shared< Grid[] >( levels_ );
   r_ = make_shared< Grid[] >( levels_ );

   for( size_t i = 0; i < nlevels; i++ )
   {
      x_[i].resize( i + 1 );
      f_[i].resize( i + 1 );
      r_[i].resize( i + 1 );
   }

   //setSmootherStencil init;
   Stencil smooth_stencil(0.25,0.25,0.25,0.25,0.25,0.0,0.0,0.0,0.0);
   setSmootherStencil( smooth_stencil );
   //setRestrictionStencil init;
   Stencil restriction_stencil(0.25,0.125,0.125,0.125,0.125,0.0625,0.0625,0.0625,0.0625);
   setRestrictionStencil( restriction_stencil );
   //setProlongationStencil init;
   Stencil interpolation_stencil(1.0,0.5,0.5,0.5,0.5,0.25,0.25,0.25,0.25);
   setProlongationStencil( interpolation_stencil );
}

MGSolver::~MGSolver()
{

}

void MGSolver::setNumberOfCycle()
{
  N_ = make_shared<size_t[11]>({4,3,3,3,4,4,4,4,4,5,6});

}

void MGSolver::setGhostLayer( size_t addy, size_t addx )
{
   for( size_t i = 0; i < levels_; ++i )
   {
      x_[i].resize( x_[i].ysize(), x_[i].xsize() + addx );
      f_[i].resize( f_[i].ysize(), f_[i].xsize() + addx );
      r_[i].resize( r_[i].ysize(), r_[i].xsize() + addx );
   }
}

void MGSolver::solveFMG(Grid& x, const Grid& f, size_t pre, size_t post)
{
        real res_cur;
  
      calcResidual(x,f,r_[levels_-1]);
      restrict(r_[levels_-1],f_[levels_-2]);
      x_[levels_ - 2].fill( 0.0 );
      for (size_t level = levels_-2;level>0;--level){
	 // Calculating the residual
         calcResidual( x_[level], f_[level], r_[level] );
         // Restriction
         restrict( r_[level], f_[level - 1] );
         // Resetting the approximate solution on the coarse grid
         x_[level - 1].fill( 0.0 );
      }
      
      // Solving the coarsest level
      smooth( x_[0], f_[0], 1 );
      //go down to coarsest level
      
      // Full multigrid solving
      for (size_t level = 1;level< levels_-1;++level){
	///  >
	/// /
	// Prolongation & error correction
	 prolongateAndCorrect( x_[level - 1], x_[level] );
	 // Post - smoothing
	 smooth( x_[level], f_[level], post );
	 
	 //V-cycle
	 solve(x_[level],f_[level],N_[level],pre,post,level+1);
	 
	 
 	 /* \
// 	     \
// 	      > */
// 	 for (size_t sublevel = level;sublevel>0;--sublevel){
// 	    // Pre - smoothing
// 	    smooth( x_[sublevel], f_[sublevel], pre );
// 	    // Calculating the residual
// 	    calcResidual( x_[sublevel], f_[sublevel], r_[sublevel] );
// 	    // Restriction
// 	    restrict( r_[sublevel], f_[sublevel - 1] );
// 	    // Resetting the approximate solution on the coarse grid
// 	    x_[sublevel - 1].fill( 0.0 );	   
// 	 }
// 	 // Solving the coarsest level
// 	 smooth(x_[0], f_[0], 1);
// 	 ///   >
// 	 ///  /
// 	 /// /
// 	 for (size_t sublevel = 1;sublevel<=level;++sublevel){
// 	    // Prolongation & error correction
// 	    prolongateAndCorrect( x_[sublevel - 1], x_[sublevel] );
// 
// 	    // Post - smoothing
// 	    smooth( x_[sublevel], f_[sublevel], post );
// 	 } 
      }
      // Prolongation & error correction on the finest grid
      prolongateAndCorrect( x_[levels_ - 2], x );

      // Post - smoothing on the finest level
      smooth( x, f, post );
      
      solve(x,f,N_[levels_-1],pre,post,levels_);

      //Printing the L2 - norm of the residual
      calcResidual( x, f, r_[levels_ - 1] );
      
       res_cur = r_[levels_ - 1].calcL2Norm();
       cout << "L2 - norm of the residual = " << setw( 12 ) << res_cur<<setw( 9 )<<endl;

  
}

void MGSolver::solve( Grid& x, const Grid& f, size_t cycles, size_t pre, size_t post ,size_t levels)
{

   real res_cur, res_old;

   for( size_t cycle = 0; cycle < cycles; ++cycle )
   {
      res_old = r_[levels - 1].calcL2Norm();

      // Pre - smoothing on the finest level
      smooth( x, f, pre );

      // Calculating the residual of the finest level
      calcResidual( x, f, r_[levels - 1] );

      // Restricting the residual
      restrict( r_[levels - 1], f_[levels - 2] );

      x_[levels - 2].fill( 0.0 );

      // Downward phase of the V- cycle
      for( size_t level = levels - 2; level > 0; --level )
      {
         // Pre - smoothing
         smooth( x_[level], f_[level], pre );

         // Calculating the residual
         calcResidual( x_[level], f_[level], r_[level] );

         // Restriction
         restrict( r_[level], f_[level - 1] );

         // Resetting the approximate solution on the coarse grid
         x_[level - 1].fill( 0.0 );
      }

      // Solving the coarsest level
      smooth( x_[0], f_[0], 1 );

      // Upward phase of the V- cycle
      for( size_t level = 1; level < levels - 1; ++level )
      {
         // Prolongation & error correction
         prolongateAndCorrect( x_[level - 1], x_[level] );

         // Post - smoothing
         smooth( x_[level], f_[level], post );
      }

      // Prolongation & error correction on the finest grid
      prolongateAndCorrect( x_[levels - 2], x );

      // Post - smoothing on the finest level
      smooth( x, f, post );
      // Printing the L2 - norm of the residual
      calcResidual( x, f, r_[levels - 1] );

      res_cur = r_[levels - 1].calcL2Norm();

      cout <<"LvL "<<setw(3)<<levels <<" Cycle " << std::setw( 3 ) << cycle + 1 << ": L2 - norm of the residual = " << setw( 12 ) << res_cur << " Convergence rate q = " << setw( 9 )
      << res_cur / res_old << "\n";
   }

}

void MGSolver::setSmootherStencil( const Stencil& stencil )
{
   A_ = stencil;
}

void MGSolver::setRestrictionStencil( const Stencil& stencil )
{
   I_ = stencil;
}

void MGSolver::setProlongationStencil( const Stencil& stencil )
{

   J_ = stencil;
}

void MGSolver::smooth( Grid& x, const Grid& f, size_t iter )
{

    //Gauss_Seidel relaxation with red-black
    #pragma omp parallel for
    for( size_t i = 1; i < x.ysize(); i++ ){
      #pragma omp parallel for
      for( size_t j = 1; j < x.xsize() - 1; j+=2 ){
	if(i!=x.ysize()-1)
	  x( i, j ) = x.hsize()*x.hsize() * A_.C * f( i, j )
		  + ( A_.N * x( i + 1, j ) + A_.S * x( i - 1, j ) + A_.W * x( i, j - 1 ) + A_.E * x( i, j + 1 ) );
	else if (j<(x.xsize()-1)/2)
	  x( i, j ) = x.hsize()*x.hsize() * A_.C * f( i, j )
		  + ( A_.N * x( i - 1, j ) + A_.S * x( i - 1, j ) + A_.W * x( i, j - 1 ) + A_.E * x( i, j + 1 ) );
	}
    }

    #pragma omp parallel for
    for( size_t i = 1; i < x.ysize(); i++ ){
      #pragma omp parallel for
      for( size_t j = 2; j < x.xsize() - 1; j+=2 ){
	if(i!=x.ysize()-1)
	  x( i, j ) = x.hsize()*x.hsize() * A_.C * f( i, j )
		  + ( A_.N * x( i + 1, j ) + A_.S * x( i - 1, j ) + A_.W * x( i, j - 1 ) + A_.E * x( i, j + 1 ) );
	else if (j<(x.xsize()-1)/2)
	  x( i, j ) = x.hsize()*x.hsize() * A_.C * f( i, j )
		  + ( A_.N * x( i - 1, j ) + A_.S * x( i - 1, j ) + A_.W * x( i, j - 1 ) + A_.E * x( i, j + 1 ) );
	}
    }

}

void MGSolver::calcResidual( const Grid& x, const Grid& f, Grid& r )
{
  #pragma omp parallel for
  for( size_t i = 1; i < x.ysize(); i++ ){
    #pragma omp parallel for
    for( size_t j = 1; j < x.xsize() - 1; j++ )
	if(i!=x.ysize()-1)
	  r( i, j ) = f( i, j ) + 1.0 / ( x.hsize()*x.hsize() ) * ( -4 * x( i, j ) + x( i + 1, j ) + x( i - 1, j ) + x( i, j + 1 ) + x( i, j - 1 ) );
	else if (j<(x.xsize()-1)/2)
	  r( i, j ) = f( i, j ) + 1.0 / ( x.hsize()*x.hsize() ) * ( -4 * x( i, j ) + x( i - 1, j ) + x( i - 1, j ) + x( i, j + 1 ) + x( i, j - 1 ) );
  }
    
}

void MGSolver::restrict( const Grid & r, Grid & f )
{
   #pragma omp parallel for
   for( size_t i = 1; i < f.ysize() ; ++i ){
      #pragma omp parallel for
      for( size_t j = 1; j < f.xsize() - 1; ++j ){
	if(i!=f.ysize()-1)
            f( i, j ) = I_.C * r( i * 2, j * 2 ) + I_.N * r( i * 2 + 1, j * 2 ) + I_.S * r( i * 2 - 1, j * 2 ) + I_.W * r( i * 2, j * 2 - 1 )
                        + I_.E * r( i * 2, j * 2 + 1 ) + I_.NW * r( i * 2 + 1, j * 2 - 1 ) + I_.NE * r( i * 2 + 1, j * 2 + 1 )
                        + I_.SW * r( i * 2 - 1, j * 2 - 1 ) + I_.SE * r( i * 2 - 1, j * 2 + 1 );
	else if (j<(f.xsize()-1)/2)
            f( i, j ) = I_.C * r( i * 2, j * 2 ) + I_.N * r( i * 2 - 1, j * 2 ) + I_.S * r( i * 2 - 1, j * 2 ) + I_.W * r( i * 2, j * 2 - 1 )
                        + I_.E * r( i * 2, j * 2 + 1 ) + I_.NW * r( i * 2 - 1, j * 2 - 1 ) + I_.NE * r( i * 2 - 1, j * 2 + 1 )
                        + I_.SW * r( i * 2 - 1, j * 2 - 1 ) + I_.SE * r( i * 2 - 1, j * 2 + 1 );

      }

   }

}

void MGSolver::prolongateAndCorrect( const Grid& e, Grid& x )
{
   #pragma omp parallel for
   for( size_t i = 1; i < e.ysize() ; ++i ){
      #pragma omp parallel for
      for( size_t j = 1; j < e.xsize() - 1; ++j ){
	if(i!=e.ysize()-1){
            x( i * 2, j * 2 ) += J_.C * e( i, j );
            x( i * 2 + 1, j * 2 ) += J_.N * e( i, j );
            x( i * 2 - 1, j * 2 ) += J_.S * e( i, j );
            x( i * 2, j * 2 - 1 ) += J_.W * e( i, j );
            x( i * 2, j * 2 + 1 ) += J_.E * e( i, j );
            x( i * 2 + 1, j * 2 - 1 ) += J_.NW * e( i, j );
            x( i * 2 + 1, j * 2 + 1 ) += J_.NE * e( i, j );
            x( i * 2 - 1, j * 2 - 1 ) += J_.SW * e( i, j );
            x( i * 2 - 1, j * 2 + 1 ) += J_.SE * e( i, j );
	}else if (j<(e.xsize()-1)/2){
	    x( i * 2, j * 2 ) += J_.C * e( i, j );
            x( i * 2 - 1, j * 2 ) += J_.S * e( i, j );
            x( i * 2, j * 2 - 1 ) += J_.W * e( i, j );
            x( i * 2, j * 2 + 1 ) += J_.E * e( i, j );
            x( i * 2 - 1, j * 2 - 1 ) += J_.SW * e( i, j );
            x( i * 2 - 1, j * 2 + 1 ) += J_.SE * e( i, j );
	}
      }
   }

}

#endif // MGSOLVER_H
