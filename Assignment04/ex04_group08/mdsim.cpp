#include <iostream>
#include "Setup.hpp"
#include <cmath>

//LINKCELL
void Linkcell(std::size_t N, std::vector<std::vector<double > >& x, std::vector<double> &min,
	std::vector<double>& max, std::vector<double> &boxlength, std::vector<std::vector<int > > &cell_num,
	std::vector<double> &cell_length, std::vector<int> &particles, std::vector<int> &cells, std::vector<int> &num_of_cell){
	
	//first check if particles are out of bound - implement periodic condition
	for(std::size_t i=0; i<N; ++i){
		for(int j=0; j<3; ++j) {
			if(x[i][j]<min[j])
				x[i][j] += boxlength[j];
			if(x[i][j]>= max[j])
				x[i][j] -= boxlength[j];
		}
	}

	//Then implement link cell
	for(std::size_t i=0; i<N; ++i){
		for(int j=0; j<3; ++j){
			cell_num[i][j] = x[i][j]/cell_length[j];
		}
		int cell_num_total = cell_num[i][0]+cell_num[i][1]*num_of_cell[0]+cell_num[i][2]*num_of_cell[1]*num_of_cell[0];
		int temp = particles[i];
		particles[i]=cells[cell_num_total];
		cells[cell_num_total] = temp;
	}
}

//FORCE CAL

void ForceCal(std::size_t N, std::vector<std::vector<double > >& x, std::vector<double> &min,
	std::vector<double>& max, std::vector<double> &boxlength, std::vector<std::vector<int > > &cell_num,
	std::vector<double> &cell_length, std::vector<int> &particles, std::vector<int> &cells, 
	std::vector<int> &num_of_cell, double r_cut, double sigma, std::vector<std::vector<double > > &F, double epsilon){


	for(std::size_t i=0; i<N; ++i){
		for(int j=0; j<3; ++j){
			//reset the force vector
			F[i][j] = 0.0;
		}
		//Treat periodic boundary condition (if neighbour box is out of domain, loop back)
		for(int k=-1;k<2;++k){
			for(int l=-1; l<2;++l){
				for(int m=-1; m<2; ++m){
					std::vector<int> cell_to_calc (3,0);//the neighbour cell we are looking at moment
					std::vector<double> offset (3,0); //this is to offset the distance between the particles later
					cell_to_calc[0] = cell_num[i][0] + k;
					cell_to_calc[1] = cell_num[i][1] + l;
					cell_to_calc[2] = cell_num[i][2] + m;
					for(int j=0; j<3; ++j){
						if(cell_to_calc[j]<0){
							cell_to_calc[j]+=num_of_cell[j];
							offset[j] = -boxlength[j];
						}
						if(cell_to_calc[j]>=num_of_cell[j]){
							cell_to_calc[j]-=num_of_cell[j];
							offset[j] = boxlength[j];
						}
					}
						int cell_to_calc_num = cell_to_calc[0]+cell_to_calc[1]*num_of_cell[0]+cell_to_calc[2]*num_of_cell[0]*num_of_cell[1];
						int k = cells[cell_to_calc_num]; //to get the particle number in the cell

						while(k!=-1) //if k=-1, then no more particles in the cell
						{
							if(k!=(int)i) //if i am not the particle in question
							{
								std::vector<double> r_ij(3,0.0);
								for(int j=0; j<3; ++j){
									r_ij[j] = x[i][j] - (x[k][j]+ offset[j]) ;
								}
								double r_ij_norm_sqr =0;
								for(int j=0; j<3; ++j){
									r_ij_norm_sqr+=pow(r_ij[j],2);
								}
								if(sqrt(r_ij_norm_sqr) <=r_cut){
									for(int j=0; j<3; ++j){
										double pow6 = pow(pow(sigma,2)/r_ij_norm_sqr,3);
										F[i][j] += 24.0 *epsilon/r_ij_norm_sqr*pow6*(2*pow6-1)*r_ij[j];
									}
								}
							}
							k = particles[k];
						}
					}
			}
		}
	}
}

void Position_update(double N, std::vector<std::vector<double > > &x, std::vector<std::vector<double > > &v, 
	std::vector<std::vector<double > > &F, double delta_t, std::vector<double> & mass){ 
	
		for(std::size_t i=0; i<N; ++i){
			for(int j=0; j<3; ++j){
				x[i][j] += delta_t*v[i][j] + F[i][j]*pow(delta_t,2)/(2*mass[i]);
			}
		}
}

void Velocity_update(double N, std::vector<std::vector<double > > &v, std::vector<std::vector<double > > &F,
	 std::vector<std::vector<double > > &F_old, double delta_t, std::vector<double> & mass){

		 for(std::size_t i=0; i<N; ++i){
			 for(int j=0; j<3; ++j){
				 v[i][j] += (F_old[i][j] + F[i][j])*delta_t/(2*mass[i]);
			 }
		 }
}

void vtk_output(int vtk_out_freq, int iteration, std::string & vtk_out_name_base, int& vtk_out_count, std::size_t N, 
		const std::vector<std::vector<double > >& x,const std::vector<std::vector<double > >& F, 
		const std::vector<std::vector<double > >& v, const std::vector<double> & mass)
{
    if (iteration % vtk_out_freq == 0){
		std::ostringstream sstm;
		std::ofstream vtk_out_stream;
		sstm << vtk_out_name_base << vtk_out_count << ".vtk";
		vtk_out_stream.open(sstm.str().c_str());

		vtk_out_stream << "# vtk DataFile Version 3.0" << std::endl;
		vtk_out_stream << "SiWiRVisFile" << std::endl;
		vtk_out_stream << "ASCII" << std::endl;
		vtk_out_stream << "DATASET UNSTRUCTURED_GRID" << std::endl;
		vtk_out_stream << "POINTS " << N << " DOUBLE" << std::endl;

		for(std::size_t i=0 ; i < N; ++i) {
			vtk_out_stream << x[i][0] <<" " << x[i][1] << " " << x[i][2] << std::endl;
		 }
		 
		vtk_out_stream << "POINT_DATA " << N<<std::endl;
		vtk_out_stream << "SCALARS mass double" << std::endl;
		vtk_out_stream << "LOOKUP_TABLE default" << std::endl;

		for(std::size_t i=0 ; i < N; ++i) {
			vtk_out_stream<< mass[i]<<std::endl;
		 }

		vtk_out_stream << "VECTORS force double" << std::endl;

		for(std::size_t i=0 ; i < N; ++i) {
			vtk_out_stream<< F[i][0] << " " << F[i][1] << " " << F[i][2] << std::endl;
		 }

		 vtk_out_stream << "VECTORS velocity double" << std::endl;

		for(std::size_t i=0 ; i < N; ++i) {
 			vtk_out_stream<< v[i][0] << " " << v[i][1] << " " << v[i][2] << std::endl;
		 }

		++(vtk_out_count);
		vtk_out_stream.close();
	}
}

int main(int argc, const char* argv[]) {
  
  	if (argc < 3)
	{
		std::cout << "No input or data file specified" << std::endl;
		return 0;
	}
	
	//Setup my prob by reading in parameters and initializing all required data
	Setup setup(argv[1], argv[2]);
	int iteration = 0;
			  		 
	vtk_output(setup.vis_space, iteration, setup.name, setup.vtk_out_count, setup.N, 
	setup.x,setup.F, setup.v, setup.mass);
		
	iteration++;

	//Linkcell algo

		//first check if particles are out of bound - implement periodic condition
		for(std::size_t i=0; i<setup.N; ++i){
			for(int j=0; j<3; ++j) {
				if(setup.x[i][j]<setup.min[j])
					setup.x[i][j] += setup.boxlength[j];
				if(setup.x[i][j]>= setup.max[j])
					setup.x[i][j] -= setup.boxlength[j];
			}
		}

		//Then implement link cell
		Linkcell(setup.N, setup.x, setup.min, setup.max, setup.boxlength, setup.cell_num, setup.cell_length, setup.particles, setup.cells, setup.num_of_cell);

		
		//Initial Force Calculation

		ForceCal(setup.N, setup.x, setup.min, setup.max, setup.boxlength, setup.cell_num, setup.cell_length, setup.particles, setup.cells, setup.num_of_cell,
			setup.r_cut, setup.sigma, setup.F, setup.epsilon);
		
		//Velocity-Verlet
		for(double t=0; t<setup.t_end; t+=setup.delta_t)
		{

			// position update
			Position_update(setup.N, setup.x, setup.v, setup.F, setup.delta_t, setup.mass);

			// reset linkcell
			fill(setup.cells.begin(), setup.cells.end(), -1);

			int j = 0;
			for (std::vector<int>::iterator it = setup.particles.begin(); it != setup.particles.end(); ++it) {
				*it = j;
				++j;
			}

			//linkcell
			Linkcell(setup.N, setup.x, setup.min, setup.max, setup.boxlength, setup.cell_num, setup.cell_length, setup.particles, setup.cells, setup.num_of_cell);

			// copy forces
			for(std::size_t i=0; i<setup.N; ++i){
				for(int j=0; j<3; ++j){
					setup.F_old[i][j] = setup.F[i][j];
				}
			}

			//Force calculation
			ForceCal(setup.N, setup.x, setup.min, setup.max, setup.boxlength, setup.cell_num, setup.cell_length, setup.particles, setup.cells, setup.num_of_cell,
			setup.r_cut, setup.sigma, setup.F, setup.epsilon);

			//update Velocites
			Velocity_update(setup.N, setup.v, setup.F, setup.F_old, setup.delta_t, setup.mass);
			
						  		 
			vtk_output(setup.vis_space, iteration, setup.name, setup.vtk_out_count, setup.N, 
			setup.x,setup.F, setup.v, setup.mass);
		
			iteration++;

		}





	return 0;

}
