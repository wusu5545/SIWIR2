#include "functions.h"
#include <sstream>
#include <math.h>

extern string name;
extern size_t vis_space;
extern size_t N;
extern vector<double> x;
extern vector<double> v;
extern boost::shared_ptr<double[]> F;
extern vector<double> mass;


void vtk_output(size_t iteration,size_t &vtk_count)
{
  if (iteration % vis_space == 0){
    ostringstream output_file;
    ofstream output;
    output_file <<"vtk/"<< name <<vtk_count<< ".vtk";
    output.open(output_file.str());
    
    output<< "# vtk DataFile Version 3.0" <<endl;
    output<< "SiWiRVisFile" <<endl;
    output<< "ASCII" <<endl;
    output<< "DATASET UNSTRUCTURED_GRID" <<endl;
    output<< "POINTS "<<N<<" DOUBLE" <<endl;
    
    for (size_t i=0;i<N;++i)
      output << x(i,0) <<" "<< x(i,1) <<" "<< x(i,2) <<endl;
        
    output<< "POINT_DATA " << N <<endl;
    output<< "SCALARS mass double" <<endl;
    output<< "LOOKUP_TABLE default" <<endl;
    
    for (size_t i=0;i<N;++i)
      output<<mass[i]<<endl;
    
    output<< "VECTORS force double" <<endl;

    for (size_t i=0;i<N;++i)
      output<<F(i,0)<<" "<<F(i,1)<<" "<<F(i,2)<<endl;
    
    output<< "VECTORS velocity double" <<endl;
    
    for (size_t i=0;i<N;++i)
      output << v(i,0)<<" "<<v(i,1)<<" "<< v(i,2)<<endl;
    
    ++vtk_count;
    output.close();
  }
}
