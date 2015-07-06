#include "functions.h"
#include "ParameterReader.h"

//parameter file
extern string name;
extern size_t vis_space;
extern double t_start;
extern double t_end;
extern double delta_t;
extern boost::shared_ptr<double[]> minxyz;
extern boost::shared_ptr<double[]> maxxyz;
extern double r_cut;
extern double epsilon;
extern double sigma;

//data file
extern vector<double> x;
extern vector<double> v;
extern boost::shared_ptr<double[]> F;
extern boost::shared_ptr<double[]> F_old;
extern boost::shared_ptr<double[]> r;
extern vector<double> mass;

//derived parameters
extern size_t N;
extern boost::shared_ptr<double[]> boxsize;
extern boost::shared_ptr<size_t[]> cells_of_xyz;
extern boost::shared_ptr<double[]> cellsize;


//link cell
extern boost::shared_ptr<int[]> cells;
extern boost::shared_ptr<size_t[]> particles;
extern boost::shared_ptr<size_t[]> cell_num;

void Initialization(string parfile,string datfile)
{
  ParameterReader reader;
  reader.readParameters(parfile);
  reader.GetParameter("name",name);
  reader.GetParameter("vis_space",vis_space);
  reader.GetParameter("t_start",t_start);
  reader.GetParameter("t_end",t_end);
  reader.GetParameter("delta_t",delta_t);
  reader.GetParameter("x_min",minxyz[0]);
  reader.GetParameter("y_min",minxyz[1]);
  reader.GetParameter("z_min",minxyz[2]);
  reader.GetParameter("x_max",maxxyz[0]);
  reader.GetParameter("y_max",maxxyz[1]);
  reader.GetParameter("z_max",maxxyz[2]);
  reader.GetParameter("r_cut",r_cut);
  reader.GetParameter("epsilon",epsilon);
  reader.GetParameter("sigma",sigma);
  
  //read dat file and initializing the data
  ifstream input(datfile);
  
  double mass_temp;
  while(input>>mass_temp)
  {
    double x_temp;
    double v_temp;
    
    mass.push_back(mass_temp);
    for (size_t i=0;i<3;++i){
      input>>x_temp;
      x.push_back(x_temp);
    }
    for (size_t i=0;i<3;++i){
      input>>v_temp;
      v.push_back(v_temp);
    }
  }
  
  //initializing all derived parameters
  N = x.size()/3;
  boxsize[0] = maxxyz[0] - minxyz[0];
  boxsize[1] = maxxyz[1] - minxyz[1];
  boxsize[2] = maxxyz[2] - minxyz[2];
  
  for (size_t i=0;i<3;++i){
    cells_of_xyz[i] = ceil(boxsize[i]/r_cut);
    cellsize[i] = boxsize[i]/cells_of_xyz[i];
  }
  
  size_t num_of_cells = cells_of_xyz[0]*cells_of_xyz[1]*cells_of_xyz[2];
  cells = boost::make_shared<int[]>(num_of_cells,-1);
  particles = boost::make_shared<size_t[]>(N);
  for (size_t i=0;i<N;++i)
    particles[i] = i;
  
  F = boost::make_shared<double[]>(N*3);
  F_old = boost::make_shared<double[]>(N*3);
  r = boost::make_shared<double[]>(N*3);
  
}
