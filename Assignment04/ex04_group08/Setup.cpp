#include "Setup.hpp" 
#include <algorithm>
  
Setup::Setup(std::string file_name, std::string input_file){
  
  //Reading in .par file and storing the data
	Parameter_Reader reader;
	reader.readParameters(file_name);
	min.resize(3);
	max.resize(3);
	reader.GetParameter("name", name);
	reader.GetParameter("vis_space", vis_space);
	reader.GetParameter("t_start", t_start);
	reader.GetParameter("t_end", t_end);
	reader.GetParameter("delta_t", delta_t);
	reader.GetParameter("x_min", min[0]);
	reader.GetParameter("y_min", min[1]);
	reader.GetParameter("z_min", min[2]);
	reader.GetParameter("x_max", max[0]);
	reader.GetParameter("y_max", max[1]);
	reader.GetParameter("z_max", max[2]);
	reader.GetParameter("r_cut", r_cut);
	reader.GetParameter("epsilon", epsilon);
	reader.GetParameter("sigma", sigma);

	//Reading in .in file and initializing the data
	std::ifstream readIn(input_file.c_str());

	double mass_temp;
	while(readIn>>mass_temp)
	{
		std::vector<double> x_temp(3,0.0);
		std::vector<double> v_temp(3,0.0);

		mass.push_back(mass_temp);
		for(int j=0; j<3; ++j) 
		{
			readIn>>x_temp[j];
		}
		for(int j=0; j<3; ++j)
		{
			readIn>>v_temp[j];
		}
		x.push_back(x_temp);
		v.push_back(v_temp);
	}

	//All derived parameters to be initialized
	N = x.size();
	boxlength.push_back(max[0] - min[0]);
	boxlength.push_back(max[1] - min[1]);
	boxlength.push_back(max[2] - min[2]);
	num_of_cell.resize(3);
	for(int i=0; i<3; ++i)
		num_of_cell[i] = ceil(boxlength[i]/r_cut);
	cell_length.resize(3);
	for(int i=0; i<3; ++i)
		cell_length[i] = boxlength[i]/num_of_cell[i];

	//For linkcell algo
	int num_cell = num_of_cell[0] * num_of_cell[1] * num_of_cell[2];
	cells.resize(num_cell);
	fill(cells.begin(), cells.end(), -1);
	
	particles.resize(N);
	int j = 0;
	for (std::vector<int>::iterator it = particles.begin(); it != particles.end(); ++it) {
		*it = j;
		++j;
	}

	cell_num.resize(N);
	for(std::size_t i=0; i<N;++i)
		cell_num[i].resize(3);

	//Resize all vectors required
	F.resize(N);
	for(std::size_t i=0; i<N; ++i)
		F[i].resize(3);
	F_old.resize(N);
	for(std::size_t i=0; i<N; ++i)
		F_old[i].resize(3);
	r_ij.resize(N);
	for(std::size_t i=0; i<N; ++i)
		r_ij[i].resize(3);

	//For output
	vtk_out_count = 0;
}

Setup::~Setup(){
}
