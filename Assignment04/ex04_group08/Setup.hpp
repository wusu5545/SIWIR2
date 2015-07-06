#ifndef Setup

#include "ParameterReader.hpp"
#include <vector>

class Setup
{   
public:
  
	//From parameter files
	std::string name;
	unsigned int vis_space;
	double t_start;
	double t_end;
	double delta_t;
	std::vector<double> min;
	std::vector<double> max;
	double r_cut;
	double epsilon;
	double sigma;
	
	//For output
	int vtk_out_count;
  
	//From input file
	std::vector<std::vector<double > > x;
	std::vector<std::vector<double > > v;
	std::vector<std::vector<double > > F;
	std::vector<std::vector<double > > F_old;
	std::vector<std::vector<double > > r_ij;
	std::vector<double> mass;

	//derived
	std::size_t N;
	std::vector<double> boxlength;
	std::vector<double> cell_length;
	std::vector<int> num_of_cell;
	//For linkcell algo
	std::vector<int> cells;
	std::vector<int> particles;
	std::vector<std::vector<int > > cell_num;

	Setup(std::string file_name, std::string input_file);
  
	~Setup();
};

#endif