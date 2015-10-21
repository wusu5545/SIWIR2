#include "functions.h"
#include <sstream>
#include <math.h>

void vtk_output(size_t vtk_step,size_t t,string &vtk_file,size_t &vtk_count,
		size_t sizex,size_t sizey,const Flags &flags,const D_Field &rho,const V_Field &v)
{
  if (t % vtk_step == 0){
    ostringstream output_file;
    ofstream output;
    output_file << vtk_file <<vtk_count<< ".vtk";
    output.open(output_file.str());
    
    output<< "# vtk DataFile Version 4.0" <<endl;
    output<< "SiwiRVisFile" <<endl;
    output<< "ASCII" <<endl;
    output<< "DATASET STRUCTURED_POINTS" <<endl;
    output<< "DIMENSIONS " <<sizex-2<<" "<<sizey-2<<" "<< 1 <<endl;
    output<< "ORIGIN 0 0 0" <<endl;
    output<< "SPACING 1 1 1" <<endl;
    output<< "POINT_DATA " << (sizex-2)*(sizey-2) <<endl<<endl;
    output<< "SCALARS flags double 1" <<endl;
    output<< "LOOKUP_TABLE default" <<endl;
    
    for (size_t j=1;j<sizey-1;++j)
      for (size_t i=1;i<sizex-1;++i)
	output<<flags(i,j)<<endl;
    output<<endl;
    
    output<< "SCALARS density double 1" <<endl;
    output<< "LOOKUP_TABLE default" <<endl;
    for (size_t j=1;j<sizey-1;++j)
      for (size_t i=1;i<sizex-1;++i)
	output<<rho(i,j)<<endl;
    output<<endl;
    
    output<< "VECTORS velocity double" <<endl;
    for (size_t j = 1;j<sizey-1;++j)
      for (size_t i = 1;i<sizex-1;++i)
	output<<(abs(v(i,j,0))<1e-08?0:v(i,j,0))<<" "<<(abs(v(i,j,1))<1e-08?0:v(i,j,1))<<" "<< 0 <<endl;
    output<<endl;
    
    ++vtk_count;
    output.close();
  }
}
