#include <math.h>
#include <omp.h>

#include "grid.h"


//Constructors
Grid::Grid(){
	// empty constructor
}

Grid::Grid(size_t l){
	size = (1<<l) + 1; // equals 2^l+1
	halfSize = size >> 1;  // equals (size-1)/2
	h_ = real(2)/(size-1);
	v_f = new float[size*(halfSize+1)];
	v_d = new double[size*(halfSize+1)];
	v_f += (size+1)*halfSize;
	v_d += (size+1)*halfSize;
	fill_f(0);
	fill_d(0);
}

//Destructor
Grid::~Grid(){
}

const Grid& Grid::operator=(const Grid& rhs){
	size = rhs.getSize();
	halfSize = size >> 1;
	h_ = rhs.hsize();
	v_f = rhs.getV_f();
	v_d = rhs.getV_d();
	return *this;
}


//Utility functions
double Grid::calcL2Norm_d() const{
	double l2norm = 0;
	for(int i = -halfSize+1; i < halfSize; ++i)
		for(int j = -halfSize+1; j < halfSize; ++j){
			real tmp = (*this)(-std::abs(i), j, 'd');
			l2norm += tmp * tmp;
		}

	l2norm /= (size)*(size);
	return sqrt(l2norm);
}

//Utility functions
float Grid::calcL2Norm_f() const{
	float l2norm = 0;
	for(int i = -halfSize+1; i < halfSize; ++i)
		for(int j = -halfSize+1; j < halfSize; ++j){
			real tmp = (*this)(-std::abs(i), j);
			l2norm += tmp * tmp;
		}

	l2norm /= (size)*(size);
	return sqrt(l2norm);
}

void Grid::fill_f(float value){
#pragma omp parallel for
	for(int i = -halfSize; i < 1; ++i){
		for(int j = -halfSize; j < halfSize+1; ++j){
			(*this)(i, j) = value;
		}
	}
}

void Grid::fill_d(double value){
#pragma omp parallel for
	for(int i = -halfSize; i < 1; ++i){
		for(int j = -halfSize; j < halfSize+1; ++j){
			(*this)(i, j, 'd') = value;
		}
	}
}

void Grid::setBoundary(){
	/*
	//  0° < phi <= 90°
	for(int i = 1; i < halfSize+1; ++i){
	(*this)(i, halfSize) = sqrt(sqrt(i*i*h_*h_ + 1)) * sin(real(0.5)*atan(i*h_));
	}
	for(int j = halfSize; j > 0; --j){
	(*this)(halfSize, j) = sqrt(sqrt(j*j*h_*h_ + 1)) * sin(real(0.5)*atan(real(1)/(j*h_)));
	}
	(*this)(halfSize, 0) = real(1)/sqrt(real(2));

	// 90° < phi <= 180°
	for(int j = -1; j > -halfSize-1; --j){
	(*this)(halfSize, j) = sqrt(sqrt(j*j*h_*h_ + 1)) * sin(real(0.5)*(M_PI+atan(real(1)/(j*h_))));
	}
	for(int i = halfSize; i > 0; --i){
	(*this)(i, -halfSize) = sqrt(sqrt(i*i*h_*h_ + 1)) * sin(real(0.5)*(M_PI+atan(-i*h_)));
	}
	*/

	(*this)(0, -halfSize) = 1.0f;
	(*this)(0, -halfSize,'d') = 1.0f;


	//180° < phi <= 270°
	for(int i = -1; i > -halfSize-1; --i){
		(*this)(i, -halfSize) = sqrt(sqrt(i*i*h_*h_ + 1)) * sin(real(0.5)*(M_PI+atan(-i*h_)));
		(*this)(i, -halfSize,'d') = sqrt(sqrt(i*i*h_*h_ + 1)) * sin(real(0.5)*(M_PI+atan(-i*h_)));
	}
	for(int j = -halfSize; j < 0; ++j){
		(*this)(-halfSize, j) = sqrt(sqrt(j*j*h_*h_ + 1)) * sin(real(0.5)*(M_PI+atan(real(-1)/(j*h_))));
		(*this)(-halfSize, j,'d') = sqrt(sqrt(j*j*h_*h_ + 1)) * sin(real(0.5)*(M_PI+atan(real(-1)/(j*h_))));
	}
	(*this)(-halfSize, 0) = real(1)/sqrt(real(2));
	(*this)(-halfSize, 0,'d') = real(1)/sqrt(real(2));

	//270° < phi <  360°
	for(int j = 1; j < halfSize+1; ++j){
		(*this)(-halfSize, j) = sqrt(sqrt(j*j*h_*h_ + 1)) * sin(real(0.5)*(2*M_PI+atan(real(-1)/(j*h_))));
		(*this)(-halfSize, j,'d') = sqrt(sqrt(j*j*h_*h_ + 1)) * sin(real(0.5)*(2*M_PI+atan(real(-1)/(j*h_))));
	}
	for(int i = -halfSize; i < 0; ++i){
		(*this)(i, halfSize) = sqrt(sqrt(i*i*h_*h_ + 1)) * sin(real(0.5)*(2*M_PI+atan(i*h_)));
		(*this)(i, halfSize,'d') = sqrt(sqrt(i*i*h_*h_ + 1)) * sin(real(0.5)*(2*M_PI+atan(i*h_)));
	}
}

void Grid::setExactSolution(){

	//(*this).setBoundary();
	(*this).fill_d(0);
	for(int j = -halfSize; j < 0; j++){	
		(*this)(0, j, 'd') = sqrt(sqrt(j*j*h_*h_));
		(*this)(0, j) = sqrt(sqrt(j*j*h_*h_));
	}

	//180° < phi <= 270°
	for(int i = -1; i > -halfSize-1; --i){
		for(int j = -halfSize; j < 0; ++j){
			(*this)(i, j, 'd') = sqrt(sqrt(i*i*h_*h_ + j*j*h_*h_)) * sin(real(0.5)*(M_PI+atan(real(i)/real(j))));
			(*this)(i, j) = sqrt(sqrt(i*i*h_*h_ + j*j*h_*h_)) * sin(real(0.5)*(M_PI+atan(real(i)/real(j))));
		}
		(*this)(i, 0, 'd') = sqrt(-i*h_)*(real(1)/sqrt(real(2)));
		(*this)(i, 0) = sqrt(-i*h_)*(real(1)/sqrt(real(2)));
	}


	//270° < phi <  360°
	for(int i = -1; i > -halfSize-1; --i){
		for(int j = 1; j < halfSize+1; ++j){
			(*this)(i, j, 'd') = sqrt(sqrt(j*j*h_*h_ + i*i*h_*h_)) * sin(real(0.5)*(real(2)*M_PI+atan(real(i)/real(j))));
			(*this)(i, j) = sqrt(sqrt(j*j*h_*h_ + i*i*h_*h_)) * sin(real(0.5)*(real(2)*M_PI+atan(real(i)/real(j))));
		}
	}

}

