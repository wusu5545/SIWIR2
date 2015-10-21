#include "assert.h"
#include "GRID.h"
#include "FUNCTIONS.h"
#include <ostream>
#include <fstream>
#include <sstream>

using namespace lbm;

void vtk_output(int vtk_out_freq, int iteration, std::string & vtk_out_name_base, int& vtk_out_count,
std::size_t xsize,std::size_t ysize, const FLAGS& flags,const D_Field& d, const V_Field& v)
{
    if (iteration % vtk_out_freq == 0){
		std::ostringstream sstm;
		std::ofstream vtk_out_stream;
		sstm << vtk_out_name_base << vtk_out_count << ".vtk";
		vtk_out_stream.open(sstm.str().c_str());

		vtk_out_stream << "# vtk DataFile Version 4.0" << std::endl;
		vtk_out_stream << "SiWiRVisFile" << std::endl;
		vtk_out_stream << "ASCII" << std::endl;
		vtk_out_stream << "DATASET STRUCTURED_POINTS" << std::endl;
		vtk_out_stream << "DIMENSIONS " << xsize-2 << " " << ysize-2 << " " << "1" << std::endl;
		vtk_out_stream << "ORIGIN " << "0 " << "0 " << "0" << std::endl;
		vtk_out_stream << "SPACING " << "1 " << "1 " << "1" << std::endl;
		vtk_out_stream << "POINT_DATA " << (xsize-2)*(ysize-2)<<std::endl<<std::endl;
		vtk_out_stream << "SCALARS flags double 1" << std::endl;
		vtk_out_stream << "LOOKUP_TABLE default" << std::endl;

	for(uint j=1 ; j < ysize-1; ++j) {
            for (uint i = 1; i < xsize-1; ++i){
                vtk_out_stream << std::fixed << flags(i,j) << std::endl;
            }
		}
		vtk_out_stream << std::endl;

		vtk_out_stream << "SCALARS density double 1" << std::endl;
		vtk_out_stream << "LOOKUP_TABLE default" << std::endl;

        for(uint j=1 ; j < ysize-1; ++j) {
            for (uint i = 1; i < xsize-1; ++i){
                vtk_out_stream << std::fixed << d(i,j) << std::endl;
            }
		}

		vtk_out_stream << std::endl;

		vtk_out_stream << "VECTORS velocity double" << std::endl;

		for (uint j = 1; j < ysize-1; ++j) {
			for (uint i = 1; i < xsize-1; ++i) {
				vtk_out_stream << std::fixed << v(i,j,0) << " " << v(i,j,1) << " " << "0"<<std::endl;
		}
		}

		++(vtk_out_count);
		vtk_out_stream.close();
	}
}

