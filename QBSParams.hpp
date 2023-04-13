struct QBSParams
	{
		//<><><><><><><><><><><><>
		// the physics parameters
		//<><><><><><><><><><><><>


		double p_central = 0.1; // central value of field p
	    double f_central = 0.14; // central value of field f


	    //potential V is,
	    // V = m_f m_f F^2 + m_p m_p P^2 + lambda_ff F^4 / 2 + lambda_fp F^2 P^2 + lambda_pp P^4 / 2

	    double m_f = 1.2; // mass term for field f
	    double m_p = 0.8; // mass term for field f
	    double q_f =-0.8; // electric coupling for field f
	    double q_p = 0.40383; // electric coupling for field f

	    // i have not tested the fourth order interactions much
	    // as long as they are kept small enough it is probably ok

	    double lambda_ff = 0.5; // ffff coefficient in potential
	    double lambda_fp =  0.3; // ffpp coefficient in potential
	    double lambda_pp = 0.6; // pppp coefficient in potential
	    double LAMBDA = 1.0; // rescales newtons constant, G_{ab} = LAMBDA^2 8 pi G T_{ab} / c^4



	    //<><><><><><><><><><><><><>
		// the numerical parameters
		//<><><><><><><><><><><><><>

	    // the amount of time spent looking for the eigenvalues
	    // this will slow down the solve a lot, try values in range 2-6 maybe?
	    int solve_iterations = 2;

	    // the overall resolution of the simulation
	    // probbaly want this at about 128000-256000 for high resolution 
	    // higher resolution is a little slower
	    int gridpoints = 32001; 

	    // this is the physical domain of the solution.
	    // if this is changed by a factor of A then the 
	    // gridsize should also be changed by a factor of A to keep the 
	    // resolution the same.
	    // can probaly try anything from tens to a few hundred
	    double outer_radius = 80.;
	};


// // Example that is very close to zero charge
// struct QBSParams
// 	{
// 		double p_central = 0.1; // central value of field p
// 	    double f_central = 0.14; // central value of field f
// 		double m_f = 1.2; // mass term for field f
// 	    double m_p = 0.8; // mass term for field f
// 	    double q_f =-0.8; // electric coupling for field f
// 	    double q_p = 0.40383; // has been finely tuned a little
// 	    double lambda_ff = 0.5; // ffff coefficient in potential
// 	    double lambda_fp =  0.3; // ffpp coefficient in potential
// 	    double lambda_pp = 0.6; // pppp coefficient in potential
// 	    int solve_iterations = 2;
// 	    int gridpoints = 32001; 
// 	    double outer_radius = 80.;
// 	};