#ifndef MGSOLVER_H
#define MGSOLVER_H

#include<boost/function.hpp>
#include<iostream>
#include<iomanip> 
#include<assert.h>
#include<boost/make_shared.hpp>
#include<boost/shared_ptr.hpp>
#include<math.h>

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
   real & operator ()( size_t i, size_t j );
   real operator ()( size_t i, size_t j ) const;

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
   CASE = 0;
}

Grid::Grid( size_t level )
{
   resize( level );
   CASE = 0;
}

Grid::Grid( const Grid& grid )
{
   ( *this ) = grid;
   CASE = 0;
}

Grid::~Grid()
{

}

// const Grid& Grid::operator=(const Grid& grid)
// {
//   return grid;
// }

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

// Grid& Grid::operator+=(const Grid& rhs)
// {
// }
// 
// Grid& Grid::operator-=(const Grid& rhs)
// {
// }

real Grid::calcL2Norm() const
{
   real res = 0;
   if( CASE == 0 )
      for( size_t i = 0; i < y_; i++ )
         for( size_t j = 0; j < x_; j++ )
         {
            res += v_[i * x_ + j] * v_[i * x_ + j];
         }
   else
      for( size_t i = 0; i < y_; i++ )
         for( size_t j = 1; j < x_ - 1; j++ )
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
   y_ = pow( 2, level ) + 1;
   v_ = make_shared< real[] >( x_ * y_ );
   h_ = real( 1.0 ) / ( pow( real( 2.0 ), level ) ); //h
}

void Grid::resize( size_t y, size_t x )
{
   x_ = x;
   y_ = y;
   CASE = 1;
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
real& Grid::operator ()( size_t i, size_t j )
{
   assert( i < y_ );
   assert( j < x_ );
   return v_[i * x_ + j];
}

real Grid::operator()( size_t i, size_t j ) const
{
   assert( i < y_ );
   assert( j < x_ );
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
   //Set Ghost layer
   void setGhostLayer( size_t addy, size_t addx );
   // Multigrid functions
   void solve( Grid & x, const Grid & f, size_t cycles, size_t pre, size_t post );
   // Utility functions
   void setSmootherStencil( const Stencil & stencil );
   void setRestrictionStencil( const Stencil & stencil );
   void setProlongationStencil( const Stencil & stencil );
private:
   // Multigrid functions
   void smooth( Grid & x, const Grid & f, size_t iter );
   void calcResidual( const Grid & x, const Grid & f, Grid & r );
   void restrict( const Grid & r, Grid & f );
   void prolongateAndCorrect( const Grid & e, Grid & x );

   // Member variables
   bool CASE;
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
   CASE = 0;
   levels_ = nlevels;
   x_ = make_shared< Grid[] >( levels_ );
   f_ = make_shared< Grid[] >( levels_ );
   r_ = make_shared< Grid[] >( levels_ );

   for( size_t i = 0; i < nlevels; i++ )
   {
      x_[i].resize( i + 1 );
      f_[i].resize( i + 1 );
      r_[i].resize( i + 1 );
   }

   Stencil stencil_temp;
   //setSmootherStencil init;
   stencil_temp.C = 0.25;
   stencil_temp.N = 0.25;
   stencil_temp.S = 0.25;
   stencil_temp.W = 0.25;
   stencil_temp.E = 0.25;
   stencil_temp.NW = 0;
   stencil_temp.NE = 0;
   stencil_temp.SW = 0;
   stencil_temp.SE = 0;
   setSmootherStencil( stencil_temp );
   //setRestrictionStencil init;
   stencil_temp.C = 0.25;
   stencil_temp.N = 0.125;
   stencil_temp.S = 0.125;
   stencil_temp.W = 0.125;
   stencil_temp.E = 0.125;
   stencil_temp.NW = 0.0625;
   stencil_temp.NE = 0.0625;
   stencil_temp.SW = 0.0625;
   stencil_temp.SE = 0.0625;
   setRestrictionStencil( stencil_temp );
   //setProlongationStencil init;
   stencil_temp.C = 1;
   stencil_temp.N = 0.5;
   stencil_temp.S = 0.5;
   stencil_temp.W = 0.5;
   stencil_temp.E = 0.5;
   stencil_temp.NW = 0.25;
   stencil_temp.NE = 0.25;
   stencil_temp.SW = 0.25;
   stencil_temp.SE = 0.25;
   setProlongationStencil( stencil_temp );
}

MGSolver::~MGSolver()
{

}

void MGSolver::setGhostLayer( size_t addy, size_t addx )
{
   CASE = 1;
   for( size_t i = 0; i < levels_; ++i )
   {
      x_[i].resize( x_[i].ysize(), x_[i].xsize() + addx );
      f_[i].resize( f_[i].ysize(), f_[i].xsize() + addx );
      r_[i].resize( r_[i].ysize(), r_[i].xsize() + addx );
   }
}

void MGSolver::solve( Grid& x, const Grid& f, size_t cycles, size_t pre, size_t post )
{

   real res_cur, res_old;

   for( size_t cycle = 0; cycle < cycles; ++cycle )
   {
      res_old = r_[levels_ - 1].calcL2Norm();

      // Pre - smoothing on the finest level
      smooth( x, f, pre );

      // Calculating the residual of the finest level
      calcResidual( x, f, r_[levels_ - 1] );

      // Restricting the residual
      restrict( r_[levels_ - 1], f_[levels_ - 2] );

      x_[levels_ - 2].fill( 0.0 );

      // Downward phase of the V- cycle
      for( size_t level = levels_ - 2; level > 0; --level )
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
      for( size_t level = 1; level < levels_ - 1; ++level )
      {
         // Prolongation & error correction
         prolongateAndCorrect( x_[level - 1], x_[level] );

         // Post - smoothing
         smooth( x_[level], f_[level], post );
      }

      // Prolongation & error correction on the finest grid
      prolongateAndCorrect( x_[levels_ - 2], x );

      // Post - smoothing on the finest level
      smooth( x, f, post );
      // Printing the L2 - norm of the residual
      calcResidual( x, f, r_[levels_ - 1] );

      res_cur = r_[levels_ - 1].calcL2Norm();

      cout << " Cycle " << std::setw( 3 ) << cycle + 1 << ": L2 - norm of the residual = " << setw( 12 ) << res_cur << " Convergence rate q = " << setw( 9 )
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

   for( size_t k = 0; k < iter; k++ )
   {
      //assure boundary conditions
      if( CASE == 1 )
      {
         for( size_t i = 1; i < x.ysize() - 1; ++i )
         {
            x( i, 0 ) = x( i, 1 ) - x.hsize();
            x( i, x.xsize() - 1 ) = x( i, x.xsize() - 2 ) - x.hsize();
         }
      }

      //Gauss_Seidel relaxation with red-black
      for( size_t i = 1; i < x.ysize() - 1; i++ )
      {
         for( size_t j = 1; j < x.xsize() - 1; j++ )
         {
            if( ( i + j ) % 2 == 0 )
            {
               x( i, j ) = x.hsize()*x.hsize() * A_.C * f( i, j )
                        + ( A_.N * x( i + 1, j ) + A_.S * x( i - 1, j ) + A_.W * x( i, j - 1 ) + A_.E * x( i, j + 1 ) );
            }
         }
      }

      for( size_t i = 1; i < x.ysize() - 1; i++ )
      {
         for( size_t j = 1; j < x.xsize() - 1; j++ )
         {
            if( ( i + j ) % 2 == 1 )
            {
               x( i, j ) = x.hsize()*x.hsize() * A_.C * f( i, j )
                        + ( A_.N * x( i + 1, j ) + A_.S * x( i - 1, j ) + A_.W * x( i, j - 1 ) + A_.E * x( i, j + 1 ) );
            }
         }
      }

      //assure boundary conditions
      if( CASE == 1 )
      {
         for( size_t i = 1; i < x.ysize() - 1; ++i )
         {
            x( i, 0 ) = x( i, 1 ) - x.hsize();
            x( i, x.xsize() - 1 ) = x( i, x.xsize() - 2 ) - x.hsize();
         }
      }
   }

}

void MGSolver::calcResidual( const Grid& x, const Grid& f, Grid& r )
{
   for( size_t i = 1; i < x.ysize() - 1; i++ )
      for( size_t j = 1; j < x.xsize() - 1; j++ )
         r( i, j ) = f( i, j ) + 1.0 / ( x.hsize()*x.hsize() ) * ( -4 * x( i, j ) + x( i + 1, j ) + x( i - 1, j ) + x( i, j + 1 ) + x( i, j - 1 ) );
}

void MGSolver::restrict( const Grid & r, Grid & f )
{

   for( size_t i = 1; i < f.ysize() - 1; ++i )
   {
      for( size_t j = 1; j < f.xsize() - 1; ++j )
      {
         if( CASE == 0 )
            f( i, j ) = I_.C * r( i * 2, j * 2 ) + I_.N * r( i * 2 + 1, j * 2 ) + I_.S * r( i * 2 - 1, j * 2 ) + I_.W * r( i * 2, j * 2 - 1 )
                        + I_.E * r( i * 2, j * 2 + 1 ) + I_.NW * r( i * 2 + 1, j * 2 - 1 ) + I_.NE * r( i * 2 + 1, j * 2 + 1 )
                        + I_.SW * r( i * 2 - 1, j * 2 - 1 ) + I_.SE * r( i * 2 - 1, j * 2 + 1 );
         else
            f( i, j ) = I_.C * r( i * 2, j * 2 - 1 ) + I_.N * r( i * 2 + 1, j * 2 - 1 ) + I_.S * r( i * 2 - 1, j * 2 - 1 ) + I_.W * r( i * 2, j * 2 - 1 - 1 )
                        + I_.E * r( i * 2, j * 2 + 1 - 1 ) + I_.NW * r( i * 2 + 1, j * 2 - 1 - 1 ) + I_.NE * r( i * 2 + 1, j * 2 + 1 - 1 )
                        + I_.SW * r( i * 2 - 1, j * 2 - 1 - 1 ) + I_.SE * r( i * 2 - 1, j * 2 + 1 - 1 );
      }

   }

}

void MGSolver::prolongateAndCorrect( const Grid& e, Grid& x )
{

   for( size_t i = 1; i < e.ysize() - 1; ++i )
   {
      for( size_t j = 1; j < e.xsize() - 1; ++j )
      {
         if( CASE == 0 )
         {
            x( i * 2, j * 2 ) += J_.C * e( i, j );
            x( i * 2 + 1, j * 2 ) += J_.N * e( i, j );
            x( i * 2 - 1, j * 2 ) += J_.S * e( i, j );
            x( i * 2, j * 2 - 1 ) += J_.W * e( i, j );
            x( i * 2, j * 2 + 1 ) += J_.E * e( i, j );
            x( i * 2 + 1, j * 2 - 1 ) += J_.NW * e( i, j );
            x( i * 2 + 1, j * 2 + 1 ) += J_.NE * e( i, j );
            x( i * 2 - 1, j * 2 - 1 ) += J_.SW * e( i, j );
            x( i * 2 - 1, j * 2 + 1 ) += J_.SE * e( i, j );
         }
         else
         {
            x( i * 2, j * 2 - 1 ) += J_.C * e( i, j );
            x( i * 2 + 1, j * 2 - 1 ) += J_.N * e( i, j );
            x( i * 2 - 1, j * 2 - 1 ) += J_.S * e( i, j );
            x( i * 2, j * 2 - 1 - 1 ) += J_.W * e( i, j );
            x( i * 2, j * 2 + 1 - 1 ) += J_.E * e( i, j );
            x( i * 2 + 1, j * 2 - 1 - 1 ) += J_.NW * e( i, j );
            x( i * 2 + 1, j * 2 + 1 - 1 ) += J_.NE * e( i, j );
            x( i * 2 - 1, j * 2 - 1 - 1 ) += J_.SW * e( i, j );
            x( i * 2 - 1, j * 2 + 1 - 1 ) += J_.SE * e( i, j );
         }
      }

   }

}

#endif // MGSOLVER_H
