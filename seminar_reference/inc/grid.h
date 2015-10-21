#ifndef GRID_H
#define GRID_H

#include <iostream>
#include <iomanip>
#include <assert.h>
#include <vector>

#include "types.h"
#include "precision.h"

class Grid
{
	public:

		//Constructors
		explicit Grid();
		explicit Grid(size_t level);
		Grid(const Grid& grid);

		//Destructor
		~Grid();

		inline void print(){
			for(int i = -halfSize; i < 0; ++i){
				for(int j = -halfSize; j < 0; ++j){
					std::cout << std::setprecision(4) << std::setw(8) << (*this)(i, j, 'd') << " ";
				}
				std::cout << std::endl;
			}
			std::cout << std::endl;
		}

		//Copy assignment
		const Grid& operator=(const Grid& grid);

		Grid& operator-=(const Grid& grid){
			for(int i = -halfSize; i <= 0; i++){
				for(int j = -halfSize; j < halfSize+1; j++){
					(*this)(i,j,'d')=(*this)(i,j,'d') - grid(i,j);
					(*this)(i,j)=(*this)(i,j) - grid(i,j);
				}
			}
			return (*this);
		};
		
		//Get functions
		inline size_t getHalfSize() const{
			return halfSize;
		};
		inline size_t getSize() const{
			return size;
		}
		inline double hsize() const{
			return h_;
		}
		inline float* getV_f() const{
			return v_f;
		}
		inline double* getV_d() const{
			return v_d;
		}

		//Access functions
		inline double& operator()(int i, int j, char d){ // with origin at (0,0)
			return v_d[i*size + j];
		};
		inline float& operator()(int i, int j){
			return v_f[i*size + j];
		}

		inline double operator()(int i, int j, char d) const{
			return v_d[i*size + j];
		};

		inline double operator()(int i, int j) const{
			return v_f[i*size + j];
		};

		//Utility functions
		double calcL2Norm_d() const;
		float calcL2Norm_f() const;
		void fill_f(float value);
		void fill_d(double value);
		void setBoundary();
		void setBoundaryNeumann();
		void resize(size_t level);
		void initF2Val(real value);
		void setExactSolution();
	private:
		//Member variables
		int halfSize;
		int size;

		//The mesh size of the grid
		double h_;

		//The values contained in the grid
		float* v_f;
		double* v_d;
};
#endif
