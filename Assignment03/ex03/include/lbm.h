#ifndef LBM_H
#define LBM_H

#include <boost/shared_ptr.hpp>
#include <boost/make_shared.hpp>

using namespace std;

namespace lbm{
  // A convenient type definition
  typedef unsigned int uint;
  
  //A convenient type definition typedef unsigned int uint;
  template <typename Type,uint Cellsize>
  class Grid{
    public:
      inline Grid():xsize_(0),ysize_(0),data_(boost::make_shared<Type[]>(0)){};
      inline Grid(uint xsize,uint ysize):xsize_(xsize),ysize_(ysize)
					,data_(boost::make_shared<Type[]>(Cellsize*xsize*ysize,0)){};
      inline ~Grid(){};
      
      inline Type& operator()(uint x,uint y,uint f)
      {
	#if use_assert == 1
	assert(x < xsize_ && y < ysize_ && f < Cellsize);
	#endif
	return data_[y*xsize_*Cellsize + x*Cellsize + f];
      }
      inline Type operator()(uint x,uint y,uint f) const
      {
	#if use_assert == 1
	assert(x < xsize_ && y < ysize_ && f < Cellsize);
	#endif
	return data_[y*xsize_*Cellsize + x*Cellsize + f];
      }
      
      inline void reset()
      {
	data_ = boost::make_shared<Type[]>(Cellsize*xsize_*ysize_,0);
      }
      
      void swap(Grid& grid)/*throw()*/{
	#if use_assert == 1
	  assert(xsize_ == grid.xsize_ && ysize_ == grid.ysize_)
	#endif
// 	std::swap(y_, grid.y_);
// 	std::swap(x_, grid.x_);
	std::swap(data_, grid.data_);
      }
      
      inline void swap(Grid &a,Grid &b)
      {
	a.swap(b);
      }
  private:
      uint xsize_;	//Number of nodes in x-dimension
      uint ysize_;	//Number of nodes in y-dimension
      
      boost::shared_ptr<Type[]> data_;	//Linearized,1-dimensional representation of the 2D data grid
  };
  
  //Partial template specialization for Cellsize = 1
  template<typename Type>
  class Grid<Type,1>{
    public:
      inline Grid():xsize_(0),ysize_(0),data_(boost::make_shared<Type[]>(0)){};
      inline Grid(uint xsize,uint ysize):xsize_(xsize),ysize_(ysize)
					,data_(boost::make_shared<Type[]>(xsize*ysize,0)){};
      inline ~Grid(){};
      
      inline Type& operator()(uint x,uint y)
      {
	#if use_assert == 1
	assert(x < xsize_ && y < ysize_);
	#endif
	return data_[y*xsize_ + x];
      }
      inline Type operator()(uint x,uint y) const
      {
	#if use_assert == 1
	assert(x < xsize_ && y < ysize_);
	#endif
	return data_[y*xsize_ + x];
      }
      inline void resize(size_t x,size_t y)
      {
	xsize_ = x;
	ysize_ = y;
	data_ = boost::make_shared<Type[]>(xsize_*ysize_,0);
      }
    private:
      uint xsize_;
      uint ysize_;
      
      boost::shared_ptr<Type[]> data_;
  };
  
  //Partial template specialization for Cellsize = 0
  //No class definition => compile time error
  template<typename Type>
  class Grid<Type,0>;
  
  //Convenient type definitions
  typedef Grid<double,9>	PDF_Field;
  typedef Grid<double,2>	V_Field;
  typedef Grid<double,1>	D_Field;
  typedef Grid<size_t,1>	Flags;
  
} //namespace lbm

#endif//LBM_H
