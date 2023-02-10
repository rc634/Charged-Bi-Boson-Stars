#include <iostream>
#include "QBSParams.hpp"
#include "QBSSolution.hpp"
#include "QBSSolution.impl.hpp"
#include "cmath"
#include <vector>

int main()
{
	QBSParams a_QBSParams;

	std::cout.precision(15);

	double time_prev, time_new;

	std::ofstream mass_phi_file;
	mass_phi_file.open("mass_vs_phi.csv");

	std::cout << "\n><><><><><><><><><><><><><><><><><><><><><><><><\n";
    std::cout << "Running Code  \n";
    std::cout << "><><><><><><><><><><><><><><><><><><><><><><><><\n\n";

	for (int i = 0; i < 1; ++i)
	{

		std::cout << " Iteration : " << i << std::endl;
		time_prev = time(NULL);
		std::cout << " Gridsize : " << a_QBSParams.gridpoints << 
		"\n L : " << a_QBSParams.outer_radius << std::endl;
		QBSSolution a_QBS;
		a_QBS.set_initialcondition_params(a_QBSParams);
		a_QBS.main();
		time_new = time(NULL);
		std::cout << "\n\n><><><><><><><><><><><><><><><><><><><><><><><><\n"
	    << "Computation time : " << time_new - time_prev << " seconds" <<
	    "\n><><><><><><><><><><><><><><><><><><><><><><><><\n\n";
		a_QBS.output_csv();
	}

	mass_phi_file.close();

	return 0;
}
