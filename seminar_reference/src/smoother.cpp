#include <iostream>
#include <cstdlib>
#include <omp.h>

#include "smoother.h"
#include "precision.h"
#include "assert.h"
#include "grid.h"
#include "config.h"


void Smoother::rb_gauss_seidel_f(Grid& a, const Grid& f, size_t iter, size_t cycle) {

	const int halfSize = a.getHalfSize();
	
	if(halfSize == 1024){
		iter = 1;
	}

	if(halfSize >= 256){
		rb_gauss_seidel_f_RBAtOnce(a,f,iter,cycle);
		return;
	}
	

	
//	if(halfSize == 1024){
//		iter = 1;
//	}

/*
	if(halfSize < 33)
		iter += 3;
	if(cycle == 3 && halfSize == 4){
		iter += 5;
	}
*/

	//Overrelaxation: regular gauss-seidel with omega = 1 - Faster convergence for omega = 1.3
	//const real omega = 1.5; //1.4 - 0.05*cycle;
	const float omega = 1.52f;//1.4365f; //1.4 - 0.05*cycle;
	const float H = a.hsize();
	const float center = (H * H * omega) / 4.0f;
	const float outer  = 1.0f / (H * H);

	//omp_set_num_threads(6);

	//static bool RED = true;

	for(size_t c = 0; c < iter; ++c) {
		//red
#pragma omp parallel
		{
			/*if(halfSize == 1024){
			  if(RED) RED = false;
			  else RED = true;
			  }
			  if(RED || halfSize != 1024){*/
			// i < 0
			//
#pragma omp for nowait schedule(static)
			for (int i = -halfSize+1; i < 0; ++i) {
				int j = -halfSize + 2 - ((-i)&1);
				for (; j < halfSize; j+=2) {
					a(i, j) = (1-omega)*a(i,j) + center * (f(i, j) +
							outer * (a(i,j+1) + a(i,j-1) + 
								a(i+1,j) + a(i-1,j))); 
				}
			}
			//i = 0	
#pragma omp for nowait schedule(static)
			for (int j = -halfSize + 2; j < 0; j+=2) {
				a(0, j) = (1-omega)*a(0,j) + center * (f(0, j) +
						outer * (a(0,j+1) + a(0,j-1) + 
							a(-1,j) + a(-1,j))); 
			}
			//}

#pragma omp barrier
			//black

			//		if(!RED || halfSize != 1024){
			// i < 0
#pragma omp for nowait schedule(static)
			for (int i = -halfSize+1; i < 0; ++i) {
				int j = -halfSize + 1 + ((-i)&1);
				for (; j < halfSize; j+=2) {
					a(i, j) = (1-omega)*a(i,j) + center * (f(i, j) +
							outer * (a(i,j+1) + a(i,j-1) + 
								a(i+1,j) + a(i-1,j))); 
				}
			}
			// i = 0
#pragma omp for nowait schedule(static)
			for (int j = -halfSize + 1; j < 0; j+=2) {
				a(0, j) = (1-omega)*a(0,j) + center * (f(0, j) +
						outer * (a(0,j+1) + a(0,j-1) + 
							a(-1,j) + a(-1,j))); 
			}
			//}
		}

	}
}




void Smoother::rb_gauss_seidel_f_RBAtOnce(Grid& a, const Grid& f, size_t iter, size_t cycle) {

	const int halfSize = a.getHalfSize();


	//Overrelaxation: regular gauss-seidel with omega = 1 - Faster convergence for omega = 1.3
	const float omega = 1.39f;//1.4365f; //1.4 - 0.05*cycle;
	const float H = a.hsize();
	const float center = (H * H * omega) / 4.0f;
	const float outer  = 1.0f / (H * H);

	//omp_set_num_threads(6);

	//static bool RED = true;

	for(size_t c = 0; c < iter; ++c) {
		//red
#pragma omp parallel
		{
			int threadNum = omp_get_thread_num();
			int numRows   = halfSize/8;
			int startRow  = -halfSize+1 + threadNum * numRows;
			if(threadNum == 7) numRows-=2;
			

			if(threadNum != 0){
				//-1st row red
				for (int j = -halfSize+2; j < halfSize; j+=2) {
					a(startRow-1, j) = (1-omega)*a(startRow-1,j) + center * (f(startRow-1, j) +
							outer * (a(startRow-1,j+1) + a(startRow-1,j-1) + 
								a(startRow-2,j) + a(startRow,j)));
				}
			} else {
				//i = 0 red	
				for (int j = -halfSize + 2; j < 0; j+=2) {
					a(0, j) = (1-omega)*a(0,j) + center * (f(0, j) +
							outer * (a(0,j+1) + a(0,j-1) + 
								a(-1,j) + a(-1,j))); 
				}
			}

			for (int j = -halfSize+1; j < halfSize; j+=2) {
				a(startRow, j) = (1-omega)*a(startRow,j) + center * (f(startRow, j) +
						outer * (a(startRow,j+1) + a(startRow,j-1) + 
							a(startRow-1,j) + a(startRow+1,j))); 
			}

#pragma omp barrier

			//red and black rows at once
			//runs till -1 for red and till -2 for black
			for (int i = startRow+1; i < startRow+numRows-1; ++i) {
				int j = -halfSize + 2 - ((-i)&1);
				for (; j < halfSize; j+=2) {
					a(i, j) = (1-omega)*a(i,j) + center * (f(i, j) +
							outer * (a(i,j+1) + a(i,j-1) + 
								a(i+1,j) + a(i-1,j))); 
					a(i-1, j) = (1-omega)*a(i-1,j) + center * (f(i-1, j) +
							outer * (a(i-1,j+1) + a(i-1,j-1) + 
								a(i,j) + a(i-2,j))); 
				}
			}

//only in case one thread finishes work before the other thread has even started working
#pragma omp barrier
			if(threadNum == 7){
				for (int i = startRow+numRows-1; i < startRow+numRows; ++i) {
					int j = -halfSize + 2 - ((-i)&1);
					for (; j < halfSize; j+=2) {
						a(i, j) = (1-omega)*a(i,j) + center * (f(i, j) +
								outer * (a(i,j+1) + a(i,j-1) + 
									a(i+1,j) + a(i-1,j))); 
						a(i-1, j) = (1-omega)*a(i-1,j) + center * (f(i-1, j) +
								outer * (a(i-1,j+1) + a(i-1,j-1) + 
									a(i,j) + a(i-2,j))); 
					}
				}

				// i = -1 red
				for (int j = -halfSize + 1; j < halfSize; j+=2) {
					a(-1, j) = (1-omega)*a(-1,j) + center * (f(-1, j) +
							outer * (a(-1,j+1) + a(-1,j-1) + 
								a(0,j) + a(-2,j))); 
				}
			} else {

				// i = -2 black
				for (int j = -halfSize + 2; j < halfSize; j+=2) {
					a(startRow+numRows-2, j) = (1-omega)*a(startRow+numRows-2,j) + center * (f(startRow+numRows-2, j) +
							outer * (a(startRow+numRows-2,j+1) + a(startRow+numRows-2,j-1) + 
								a(startRow+numRows-3,j) + a(startRow+numRows-1,j))); 
				}
			
			}

			// i = -1 black (bzw -2 for 7th thread)
			for (int j = -halfSize + 1; j < halfSize; j+=2) {
				a(startRow+numRows-1, j) = (1-omega)*a(startRow+numRows-1,j) + center * (f(startRow+numRows-1, j) +
						outer * (a(startRow+numRows-1,j+1) + a(startRow+numRows-1,j-1) + 
							a(startRow+numRows-2,j) + a(startRow+numRows,j)));
			}

			
			
			if(threadNum == 7){
			
				// i = -1 black
				for (int j = -halfSize + 2; j < halfSize; j+=2) {
					a(-1, j) = (1-omega)*a(-1,j) + center * (f(-1, j) +
							outer * (a(-1,j+1) + a(-1,j-1) + 
								a(-2,j) + a(0,j))); 
				}

				// i = 0 black
				for (int j = -halfSize + 1; j < 0; j+=2) {
					a(0, j) = (1-omega)*a(0,j) + center * (f(0, j) +
							outer * (a(0,j+1) + a(0,j-1) + 
								a(-1,j) + a(-1,j))); 
				}
			}

#pragma omp barrier

		}
	}
}



