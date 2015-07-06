#include "functions.h"

//parameter file
string name;
size_t vis_space;
double t_start;
double t_end;
double delta_t;
boost::shared_ptr<double[]> minxyz = boost::make_shared<double[]>(3);
boost::shared_ptr<double[]> maxxyz = boost::make_shared<double[]>(3);
double r_cut;
double epsilon;
double sigma;

//data file
vector<double> x;
vector<double> v;
boost::shared_ptr<double[]> F;
boost::shared_ptr<double[]> F_old;
boost::shared_ptr<double[]> r;
vector<double> mass;

//derived parameters
size_t N;
boost::shared_ptr<double[]> boxsize = boost::make_shared<double[]>(3);
boost::shared_ptr<size_t[]> cells_of_xyz = boost::make_shared<size_t[]>(3);
boost::shared_ptr<double[]> cellsize = boost::make_shared<double[]>(3);


//link cell
boost::shared_ptr<int[]> cells;
boost::shared_ptr<size_t[]> particles;
boost::shared_ptr<size_t[]> cell_num;

int main(int argc, char* argv[])
{
  if (argc <= 1)
  {
    cout<< "Please input .par and .dat files" <<endl;
    exit(0);
  }
  
  Initialization(argv[1],argv[2]);
  
  return 0;
}
