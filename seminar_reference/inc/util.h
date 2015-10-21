#ifndef UTIL_H
#define UTIL_H

#include <iostream>
#include <fstream>
#include <cstdio>
#include <cstring>
#include <stdlib.h>
#include <string>
#include <assert.h>

#include "grid.h"
#include "config.h"


#if PRINTMODE == 1
static void printForPlot2(const Grid& g, const std::string filename) {
	std::fstream o(filename, std::fstream::out | std::fstream::trunc);
	assert(o && "fatal error: cannot open file");
	const int halfSize = g.getHalfSize();
	const real H = g.hsize();

	//o << "# x y u(x,y)" << std::endl;
	
	real yPos = -1;
	for (int y = -halfSize; y < halfSize+1; ++y) {
		real xPos = -1;
		for (int x = -halfSize; x < halfSize+1; ++x) {
			o << xPos << " " << yPos << " " << g(-std::abs(y), x) << std::endl;
			xPos += H;
		}
		o << std::endl;
		yPos += H;
	}

	o.close();
}
#elif PRINTMODE==2

static void printForTool2(const Grid& g, const std::string filename) {
	std::fstream o(filename, std::fstream::out | std::fstream::trunc);
	assert(o && "fatal error: cannot open file");
	const int halfSize = g.getHalfSize();
	const real H = g.hsize();

	//o << "# x y u(x,y)" << std::endl;
	
	real yPos = -1;
	for (int y = -halfSize; y < halfSize+1; ++y) {
		real xPos = -1;
		for (int x = -halfSize; x < halfSize+1; ++x) {
			o << xPos << " " << yPos << " " << g(-std::abs(y), x) << std::endl;
			xPos += H;
		}
		yPos += H;
	}

	o.close();
}
#endif

#endif

