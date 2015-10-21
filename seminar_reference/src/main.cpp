#include <iostream>
#include <fstream>
#include <cstdio>
#include <cstring>
#include <stdlib.h>
#include <sys/time.h>
#include <math.h>

#include "mgsolver.h"
#include "grid.h"
#include "util.h"
#include "types.h"
#include "config.h"

void warmupCPU();

int main (int argc, const char* argv[]) {
	if (argc != 2){
		std::cerr << "Usage: " << argv[0] << " <l>" << std::endl;
		exit(EXIT_FAILURE);
	}
	
	int l = atoi(argv[1]);
	if(l < 0){
		std::cerr << "please provide a positive number!" << std::endl;
		exit(EXIT_FAILURE);
	}
	
	std::cout << "calculating result... " << std::endl;
	
	MGSolver solver;

	Grids x_ = Grids((size_t)l);
	Grids f_ = Grids((size_t)l);
	Grids r_ = Grids((size_t)l);

	Grid d = Grid(COARSEST+1);
	Grid z = Grid(COARSEST+1);

	for(int i = 1; i < l; ++i){ 
		x_[i] = Grid(i+1);
		f_[i] = Grid(i+1);
		r_[i] = Grid(i+1);
	}

	x_[l-1].setBoundary();
#if PRINTMODE == 2
		printForTool2(x_[l-1], "init.dat");
#endif

	//siwir::Timer timer;
	std::cout<<"Your Alias: "<< "//TODO:insertAlias" <<std::endl;
	struct timeval t0, t;
	warmupCPU();
	gettimeofday(&t0, NULL);
	solver.solve(x_, f_, r_, d, z, (size_t)l, PRE, POST);
	gettimeofday(&t, NULL);
	std::cout << "Wall clock time of MG execution: " <<
		((int64_t)(t.tv_sec - t0.tv_sec) * (int64_t)1000000 +
		 (int64_t)t.tv_usec - (int64_t)t0.tv_usec) * 1e-3
		<< " ms" << std::endl;
	//std::cout << "Solving took " << timer.elapsed() << " seconds!" << std::endl;
#if CALCERROR == 1
		Grid exactSolution = Grid(l);
		exactSolution.setExactSolution();
		exactSolution -= x_[l-1];
		std::cout << "error: " << std::setprecision(11) <<  exactSolution.calcL2Norm_d() << std::endl;
#endif

#if PRINTMODE == 1
		printForPlot2(x_[l-1], "solution.dat");
#elif PRINTMODE == 2
		printForTool2(x_[l-1], "solution.dat");
#endif
}

void warmupCPU(){
	bool isPrime = true;
	long primeSum = 0;
	for(int prim = 2; prim < 100000; ++prim){
		for(int test = (int)sqrt(prim); test > 1; --test){
			if(prim % test == 0){
				isPrime = false;
			}
		}
		if(isPrime)
			primeSum += prim;
		isPrime = true;
	}

	std::cout << "Fun fact: The sum over all prime numbers from 2 to 100000 is: " << primeSum << std::endl;
}

