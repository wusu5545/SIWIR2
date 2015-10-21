#ifndef GRID
#define GRID

#include <iostream>
#include <cmath>
#include <vector>
#include "assert.h"
#include <algorithm>
#include <iomanip>

using std::ostream;
using std::istream;
using std::vector;

namespace lbm {

	typedef unsigned int uint;

	template< typename Type, uint Cellsize >
	class Grid {

	friend std::ostream &operator<< (std::ostream &os, const Grid &grid)
	{
	      for(uint i = 0; i<grid.ysize_; ++i){
		for(uint j = 0; j<grid.xsize_; ++j) {
		  for(uint k = 0; k<Cellsize; ++k){
		      os <<std::setw(10)<<std::setprecision(4)<< grid(j, grid.ysize_- i -1, k) << " ";
		    }
		  os<< " ";
		}
	      os <<std::endl;
	      }
	    return os;
	}

	private:
	  uint xsize_, ysize_;
	  Type* data_;

	public:
	  //CONSTRUCTOR
	  inline Grid(uint xsize = 0, uint ysize =0) : xsize_(xsize), ysize_(ysize), data_(new Type[Cellsize*xsize*ysize]) {};

	  //COPY CONSTRUCTOR
	  inline Grid(const Grid &grid): xsize_(grid.xsize_), ysize_(grid.ysize_), data_(grid.data_) {};

	  //ASSIGNMENT OPERATORS
	  inline Grid& operator=(const Grid &grid)
	  {
	    xsize_ = grid.xsize_;
	    ysize_ = grid.ysize_;
	    data_ = grid.data_;
	  }

	  //DESTRUCTORS
	  inline ~Grid()
	  {
	    delete data_;
	  }

	  //FUNCTION CALL OPERATOR
	  inline Type& operator() (uint x, uint y, uint f)
	  {
	    assert(x < xsize_ && y< ysize_ && f<Cellsize);
	    return data_[y*xsize_*Cellsize + x*Cellsize + f];
	  }

	  inline Type operator() (uint x, uint y, uint f) const
	  {
	     assert(x < xsize_ && y< ysize_ && f<Cellsize);
	    return data_[y*xsize_*Cellsize + x*Cellsize + f];
	  }

	  //SWAP MEMBER FUNCTION
	  void swap(Grid& grid)
	  {
	    assert(grid.xsize_ == xsize_ && grid.ysize_ == ysize_);
	    std::swap(xsize_, grid.xsize_);
	    std::swap(ysize_, grid.ysize_);
	    std::swap(data_, grid.data_);
	  }

	  //RESET
	  void reset()
	  {
	    std::fill(data_, data_+ysize_*xsize_*Cellsize, 0);
	  }

	};

	//CLASS PARTIAL TEMPLATE SPECIALIZATION for Cellsize = 1
	template <typename Type>
	class Grid<Type,1>
	{

	  friend std::ostream &operator<< (std::ostream &os, const Grid &grid)
	{
	      for(uint i = 0; i<grid.ysize_; ++i){
		for(uint j = 0; j<grid.xsize_; ++j) {
		      os << grid(j, grid.ysize_ - i -1) << " ";
		}
	      os <<std::endl;
	      }
	    return os;

	}

	  private:
	  uint xsize_, ysize_;
	  Type* data_; //an array containing my 9 data pts corresponding to my 9 directions from a cell.

	  public:
	  //CONSTRUCTOR
	  inline Grid(uint xsize = 0, uint ysize =0) : xsize_(xsize), ysize_(ysize), data_(new Type[xsize*ysize]) {};

	  //COPY CONSTRUCTOR
	  inline Grid(const Grid &grid): xsize_(grid.xsize_), ysize_(grid.ysize_), data_(grid.data_) {};

	  //ASSIGNMENT OPERATORS
	  inline Grid& operator=(const Grid &grid)
	  {
	    xsize_ = grid.xsize_;
	    ysize_ = grid.ysize_;
	    data_ = grid.data_;
	  }

	  //DESTRUCTORS
	  inline ~Grid()
	  {
	    delete data_;
	  }

	  //FUNCTION CALL OPERATOR
	  inline Type& operator() (uint x, uint y)
	  {
	    assert(x < xsize_ && y< ysize_);
	    return data_[y*xsize_ + x];
	  }

	  inline Type operator() (uint x, uint y) const
	  {
	    assert(x < xsize_ && y< ysize_);
	    return data_[y*xsize_ + x];
	  }

	};

	//CLASS PARTIAL TEMPLATE SPECIALIZATION for Cellsize = 0
	template <typename Type>
	class Grid<Type,0>;

	typedef Grid<double,9> PDF_Field;
	typedef Grid<double,2> V_Field;
	typedef Grid<double,1> D_Field;
	typedef Grid<uint,1> FLAGS;
}


#endif