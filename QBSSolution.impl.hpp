/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#if !defined(QBSSOLUTION_HPP_)
#error "This file should only be included through QBSSolution.hpp"
#endif

#ifndef QBSSOLUTION_IMPL_HPP_
#define QBSSOLUTION_IMPL_HPP_

#include <fstream>

QBSSolution::QBSSolution()
{
}

void QBSSolution::main()
{
    // lets goooooo

    w_IB_bistar(repeats);
    rk4(WP,WF);
    fix();
    mid_int = find_midint();
    force_flat(mid_int);

    //std::cout << "MIDINT : " << mid_int << ", MID_R : " << dx*mid_int << std::endl;

    rk4_asymp(WP,WF,mid_int);

    OM_INF = omega[gridsize-1];
    PSI_INF = psi[gridsize-1];
    Z_INF = z[gridsize-1];
    //std::cout << "Omega inf = " << OM_INF << ",\nPsi inf = " << PSI_INF << "\nZ inf = " << Z_INF << std::endl;
    
    for (int iteration=0; iteration < 2; iteration++)
    {
        //std::cout << "#####################\n> iter : " << iteration << std::endl;
        middle_w /= OM_INF;
        OMC /= OM_INF;
        ZC = (ZC-Z_INF)/OM_INF; 
        PSC /= PSI_INF;

        w_IB_bistar(repeats+iteration);
        rk4(WP,WF);
        fix();
        mid_int = find_midint();
        //std::cout << "LAST MIDR " << mid_int*dx << std::endl;  
        force_flat(mid_int);
        rk4_asymp(WP,WF,mid_int);
        OM_INF = omega[gridsize-1];
        PSI_INF = psi[gridsize-1];
        Z_INF = z[gridsize-1];
        //std::cout << "Omega inf = " << OM_INF << ",\nPsi inf = " << PSI_INF << "\nZ inf = " << Z_INF << std::endl;
    }

    std::cout << "\n><><><><><><><><><><><><><><><><><><><><><><><><\n";
    std::cout << "Asymptotic Values \n";
    std::cout << "><><><><><><><><><><><><><><><><><><><><><><><><\n";

    std::cout << "\n Omega inf = " << OM_INF << 
    "\n Psi inf = " << PSI_INF << 
    "\n Z inf = " << Z_INF <<
    "\n w_p = " << WP <<
    "\n w_f = " << WF 
    << std::endl;

    middle_w /= OM_INF;
    OMC /= OM_INF;
    ZC = (ZC-Z_INF)/OM_INF; 
    PSC /= PSI_INF;

    w_IB_bistar(3*repeats);
    rk4(WP,WF);
    fix();
    mid_int = find_midint();
    force_flat(mid_int);
    rk4_VAC(WP,WF,mid_int);
    

  	// lets make the diagnostic variables

    double rad, rawrad, ricc, rho_density, ddpsi, Ttt, Noether=0., Leccy_charge=0., NP, NF, dipole=0., HAM=0.;
    double NPtot=0., NFtot=0.;
    for (int j=0; j<gridsize; j++)
    {        
        rad = sqrt(pow(radius_array[j],2)+10e-14);
        rawrad = radius_array[j];
        ddpsi = DPSI_RHS(rawrad,p[j],f[j],z[j],dp[j],df[j],dz[j],psi[j],dpsi[j],omega[j],WP,WF);
        Ttt = p[j]*p[j]*pow((WP + q_p*z[j]),2) + f[j]*f[j]*pow((WF + q_f*z[j]),2);
        Ttt += omega[j]*omega[j]*(dp[j]*dp[j] +df[j]*df[j])/psi[j]/psi[j];
        Ttt += omega[j]*omega[j]*V(p[j],f[j]) + 0.5*dz[j]*dz[j]/psi[j]/psi[j];
        ricc = (2.*rawrad*dpsi[j]*dpsi[j]-4.*psi[j]*(2.*dpsi[j]+ ddpsi*rawrad)  )/(rad*(pow(psi[j],4)));
        rho_density = Ttt*pow(omega[j],-2);
        NP = (-8.*M_PI*p[j]*p[j]*(WP + q_p*z[j])*rawrad*rawrad*pow(psi[j],3)*omega[j]);
        NF = (-8.*M_PI*f[j]*f[j]*(WF + q_f*z[j])*rawrad*rawrad*pow(psi[j],3)*omega[j]);
        NPtot += NP*dx;
        NFtot += NF*dx;
        N[j] = q_p*NP + q_f*NF;
        ham[j] = (16.*M_PI*rho_density - ricc)/(fabs(ricc) + 10e-14);
        HAM += ham[j]*ham[j]*dx*rawrad/rad;
        radial_pressure[j] = psi[j]*psi[j]*(pow(p[j]*(WP+q_p*z[j]),2) + pow(f[j]*(WF+q_f*z[j]),2))/omega[j]/omega[j];
        radial_pressure[j] += df[j]*df[j] + dp[j]*dp[j] - V(f[j],p[j])*psi[j]*psi[j] - 0.5*dz[j]*dz[j]/omega[j]/omega[j];
        radial_pressure[j] *= 4.*M_PI*rawrad*rawrad*pow(psi[j],3); // the volume element on the surface of a sphere
        Noether += (NP+NF)*dx;
        Leccy_charge += (NP*q_p + NF*q_f)*dx;
        dipole += rad*(NP*q_p + NF*q_f)*dx;
    }

    ADM_MASS = -radius_array[gridsize-100.]*radius_array[gridsize-100.]*psi[gridsize-100.]*dpsi[gridsize-100.];

    
    std::cout << " AMD Mass : " << ADM_MASS ;
    std::cout << "\n\n><><><><><><><><><><><><><><><><><><><><><><><><\n";
    std::cout << "Charges \n";
    std::cout << "><><><><><><><><><><><><><><><><><><><><><><><><\n" << 
    "\n Noether Charge of p : " << NPtot << 
    "\n Noether Charge of f : " << NFtot << 
    "\n Total Noether Charge : " << Noether << 
    "\n Electric Charge of p: " << NPtot*q_p << 
    "\n Electric Charge of f: " << NFtot*q_f << 
    "\n Total Electric Charge : " << Leccy_charge << 
    "\n Charge Cancellation Factor  : " << 1.-fabs(Leccy_charge)/fabs(NPtot*q_p-NFtot*q_f) << 
    "\n Qtotal/Ntotal : " << Leccy_charge/Noether <<
    "\n Electric 'radial dipole' : " << dipole << std::endl;


    std::cout << "\n\n><><><><><><><><><><><><><><><><><><><><><><><><" <<
    "\nIntegrated Normalised Ham - Error Measure" << 
    "\n><><><><><><><><><><><><><><><><><><><><><><><><\n" <<
    "\n H2 Norm  : " << HAM << "\n";

}

// initiaalises the 5 filed variables with their central values
void QBSSolution::initialise()
{
  	p[0] = PC;
    f[0] = FC;
    z[0] = ZC;
  	omega[0] = OMC;
  	psi[0] = PSC;
  	// field gradients must be zero at zero radius for physical solutions
    dp[0] = 0.;
    df[0] = 0.;
  	dz[0] = 0.;
  	dpsi[0] = 0.;
    radius_array[0] = 0.;
}

// if the scalar field diverges it rounds down to a value wiht the same sign. This is to not affect axis crossings function (crossings) and let it accurately deal with inf/nan
void QBSSolution::fix()
{
    int bork_i=gridsize-1;
  	bool borked = false; // turns true if function gets over twice as large as central (r=0) value
  	double ptruncation, ztruncation, ftruncation;
  	for (int i = 0; i < gridsize; ++i)
  	{
    		if (not borked)
    		{
      			// if (fabs(p[i])> 1.1*p[0])
      			// {
        		// 		borked = true;
        		// 		truncation = p[i];
                //         ftruncation = f[i];
                //         ztruncation = z[i];
                //         fbork_minus_zbork = p[i]-z[i];
      			// }
                if (std::isnan(p[i]))
                {
                    borked = true;
                    bork_i = i*7/8;
                    ptruncation = p[bork_i];
                    ftruncation = f[bork_i];
                    ztruncation = z[bork_i];
                }
    		}
    		
  	}

    if (borked)
    {
        for (int i = bork_i; i < gridsize; ++i)
        {
            p[i] = ptruncation;    
            f[i] = ftruncation;
            z[i] = ztruncation;
        }
    }

   
}

// sets scalar field and gradient to zero after the point decided by function find_midint
void QBSSolution::force_flat(const int iter_crit)
{
  	for (int i = iter_crit; i < gridsize; ++i)
  	{
    		p[i] = 0.;
    		f[i] = 0.;
            dp[i] = 0.;
            df[i] = 0.;
  	}
}

// counts how many times the function crosses the axis
int QBSSolution::pcrossings()
{
    int number=0;
    for (int i = 0; i < gridsize-1; ++i)
    {
            if (p[i]*p[i+1]<0)
            {
                   number += 1;
            }
    }
    return number;
}

int QBSSolution::fcrossings()
{
    int number=0;
    for (int i = 0; i < gridsize-1; ++i)
    {
            if (f[i]*f[i+1]<0)
            {
                   number += 1;
            }
    }
    return number;
}

void QBSSolution::make_vacuum_outside_star(const double w_)
{
    mid_int = find_midint(); // the gridpoint to cutoff f
    //std::cout << "f_turning_point_index : " << mid_int << std::endl;
    //std::cout << "f_turning_point_radius : " << radius_array[mid_int] << std::endl;

    force_flat(mid_int);

    //rk4_VAC(w_,mid_int+2); // start early to fit a nice gradient fit
}

// finds the integer index at which the scalar field needs to be truncated
// int QBSSolution::find_midint()
// {
//   	int crossings = 0, mid_int;

//   	// follow the correct amount of crossings
//   	for (int i = 0; i < gridsize-1; ++i)
//   	{
//     		if (crossings == EIGEN)
//     		{
//       			mid_int = i;
//       			break;
//     		}
//     		if (f[i]*f[i+1] < 0.)
//     		{
//     			   crossings += 1;
//     		}
//   	}

//   	// climb the hill if crossings != 0, this will exit loop immediately if there is no crossings
//   	for (int i = mid_int+5; i < gridsize-1; ++i)
//   	{
//     		if (fabs(f[i+1]) < fabs(f[i]))
//     		{
//       			mid_int = i;
//       			break;
//     		}
//   	}

//   	for (int i = mid_int+5; i < gridsize-1; ++i)
//   	{
//     		if (fabs(f[i+1]) > fabs(f[i]))
//     		{
//       			mid_int = i;
//       			//std::cout << "Truncation error: " << f[i]/PC << std::endl;
//       			return mid_int;
//     		}
//   	}
//   	return gridsize-1;
// }

int QBSSolution::find_midint()
{   
    for (int i = 200; i <gridsize-2; i++)
    {
        if (p[i+1]*p[i+1]+f[i+1]*f[i+1] > p[i-1]*p[i-1]+f[i-1]*f[i-1]) return i;
    }   
    return gridsize-1;
}

// int QBSSolution::find_midint()
// {
//     int mid_int=-1;
//     double smallest_ffpluszz=p[0]*p[0]+z[0]*z[0];
//     fix();
    
//     for (int i = 200; i <gridsize-1; i++)
//     {
//         if (p[i]*p[i]+z[i]*z[i] < smallest_ffpluszz)
//             {
//                 if (isnan(p[i]) or isnan(z[i]))
//                 {
//                     break;
//                 }
//                 smallest_ffpluszz = p[i]*p[i]+z[i]*z[i];
//                 mid_int=i;
//             }
//     }
//     return mid_int;
// }

// finds an eigenvalue with a lot of nodes, (20 + desired eigenstate) by defualt
double QBSSolution::find_W()
{
  	int eigenstate;
  	double W_ = 1.;
  	while (true)
  	{
        //std::cout << " w : " << W_ << " , ";
    		initialise();
    		rk4(W_,-W_);
    		fix();
    		eigenstate = pcrossings();
        //std::cout << "eigenstate : " << eigenstate << std::endl;
    		if (eigenstate >= 4 + EIGEN){ return W_;}
    		W_*=2.;

  	}
}

double QBSSolution::find_W_soliton()
{
    bool crossed;
  	double W_ = 1.;
  	while (true)
  	{
    		initialise();
    		rk4(W_,-W_);
    		fix();
    		crossed = soliton_eigen();
    		if (crossed) return W_;
    		W_*=2.;
  	}
}

// calculates the lower limit for eigenvalue to be used in interval bisection
double QBSSolution::w_min(const double W_)
{
  	int accuracy = 400, eigenstate;
  	double w_, lower_w_;
  	for (int i=0; i<accuracy; i++)
  	{
    		w_ = W_*((double)(accuracy - i))/((double)accuracy);
    		initialise();
    		rk4(w_,-w_);
    		fix();
        //std::cout << "lower w : " << w_ << ", " ;
    		eigenstate = pcrossings();
        //std::cout << "eigenstate : " << eigenstate << std::endl;
    		if (eigenstate <= EIGEN)
    		{
                //std::cout << "Found Lower Bound on w !" << std::endl;
      			lower_w_ = w_;
      			return lower_w_;
    		}
  	}
  	return 0.;
}


// calculates the upper value for eigenvalue to be used in interval bisection
double QBSSolution::w_max(const double W_, const double lower_w_)
{
  	int accuracy = 400, eigenstate;
  	double w_, upper_w_;
  	for (int i=0; i<accuracy; i++)
  	{
    		w_ = lower_w + W_*((double)i)/((double)accuracy);
    		initialise();
    		rk4(w_,-w_);
    		fix();
        //std::cout << "max w : " << w_ << ", ";
    		eigenstate = pcrossings();
        //std::cout << "eigenstate : " << eigenstate << std::endl;
    		if (eigenstate > EIGEN )
    		{
                //std::cout << "Found Upper Bound on w !" << std::endl;
      			upper_w_ = w_;
      			return upper_w_;
    		}
  	}
  	return 0.;
}


// takes in an upper and lower eigenvalue and uses interval bisection to find the solution inbetween
double QBSSolution::w_IB(double lower_w_, double upper_w_)
{
  	int iter = 0, itermax, eigenstate, decimal_places_of_omega = 25.;
  	double middle_w_;

  	itermax = (int)((log(upper_w_)+decimal_places_of_omega*log(10.))/log(2.)); //calculate number of bisections needed (simple pen and paper calculation)
  	while (true)
  	{
    		iter++;
    		middle_w_ = 0.5*(upper_w_+lower_w_);
    		initialise();
    		rk4(middle_w_,-middle_w_);
    		fix();
        //std::cout << "IB w : " << middle_w_ << ", " ;
    		eigenstate = pcrossings();
        //std::cout << "eigenstate : " << eigenstate << std::endl;
    		if (eigenstate > EIGEN)
    		{
    			   upper_w_ = middle_w_;
    		}
    		else
    		{
    			   lower_w_ = middle_w_;
    		}
    		if ((upper_w_-lower_w_)<w_tolerance) return upper_w_;
        if (upper_w_==lower_w_) return upper_w_;
        if (iter>65) return upper_w_;
  	}
}


// takes in an upper and lower eigenvalue and uses interval bisection to find the solution inbetween
double QBSSolution::w_IB_bistar2()
{
    int iter = 0, itermax, f_num_cross, p_num_cross, decimal_places_of_omega = 25.;
    double middle_w_p, upper_w_p = 10., lower_w_p = 0.;
    double middle_w_f, upper_w_f = 10., lower_w_f = 0.;


    itermax = (int)((log(upper_w_p)+decimal_places_of_omega*log(10.))/log(2.)); //calculate number of bisections needed (simple pen and paper calculation)
    while (true)
    {
            iter++;
            middle_w_p = 0.5*(upper_w_p+lower_w_p);
            middle_w_f = 0.5*(upper_w_f+lower_w_f);
            initialise();
            rk4(middle_w_p,middle_w_f);
            fix();
            p_num_cross = pcrossings();
            f_num_cross = fcrossings();
            if (p_num_cross > 0)
            {
                   upper_w_p = middle_w_p;
            }
            else
            {
                   lower_w_p = middle_w_p;
            }

            initialise();
            rk4(middle_w_p,middle_w_f);
            fix();
            p_num_cross = pcrossings();
            f_num_cross = fcrossings();

            if (f_num_cross > 0)
            {
                   upper_w_f = middle_w_f;
            }
            else
            {
                   lower_w_f = middle_w_f;
            }
        //std::cout << "wp and wf : " << middle_w_p << " & "<< middle_w_f << std::endl;

        if (iter>65) break;
    }



    WP = middle_w_p;
    WF = middle_w_f;
}



// takes in an upper and lower eigenvalue and uses interval bisection to find the solution inbetween
double QBSSolution::w_IB_bistar(const double repeats_param)
{
    int iter = 0, itermax, f_num_cross, p_num_cross, decimal_places_of_omega = 25.;
    double middle_w_p, upper_w_p = 10., lower_w_p = 0.;
    double middle_w_f, upper_w_f = 10., lower_w_f = 0.;


    itermax = (int)((log(upper_w_p)+decimal_places_of_omega*log(10.))/log(2.)); //calculate number of bisections needed (simple pen and paper calculation)
    while (true)
    {
            iter++;
            middle_w_p = 0.5*(upper_w_p+lower_w_p);
            middle_w_f = 0.5*(upper_w_f+lower_w_f);
            initialise();
            rk4(middle_w_p,middle_w_f);
            fix();
            p_num_cross = pcrossings();
            f_num_cross = fcrossings();
            if (p_num_cross > 0)
            {
                   upper_w_p = middle_w_p;
            }
            else
            {
                   lower_w_p = middle_w_p;
            }

            // if (f_num_cross > 0)
            // {
            //        upper_w_f = middle_w_f;
            // }
            // else
            // {
            //        lower_w_f = middle_w_f;
            // }
        // std::cout << "wp and wf : " << middle_w_p << " & "<< middle_w_f << std::endl;

        if (iter>65) break;
    }

    // std::cout << "#######################" << std::endl;
    // std::cout << "swapping over" << std::endl;
    // std::cout << "#######################" << std::endl;

    iter=0;
    while (true)
    {
            iter++;
            middle_w_p = 0.5*(upper_w_p+lower_w_p);
            middle_w_f = 0.5*(upper_w_f+lower_w_f);
            initialise();
            rk4(middle_w_p,middle_w_f);
            fix();
            p_num_cross = pcrossings();
            f_num_cross = fcrossings();
            // if (p_num_cross > 0)
            // {
            //        upper_w_p = middle_w_p;
            // }
            // else
            // {
            //        lower_w_p = middle_w_p;
            // }

            if (f_num_cross > 0)
            {
                   upper_w_f = middle_w_f;
            }
            else
            {
                   lower_w_f = middle_w_f;
            }
        // std::cout << "wp and wf : " << middle_w_p << " & "<< middle_w_f << std::endl;

        if (iter>65) break;
    }

    // std::cout << "#######################" << std::endl;
    // std::cout << "swapping over" << std::endl;
    // std::cout << "#######################" << std::endl;

    for (int q=0; q<repeats_param; q++)
    {
        iter=0;
        lower_w_p = middle_w_p*0.5;
        upper_w_p = middle_w_p*2.;
        while (true)
        {
                iter++;
                middle_w_p = 0.5*(upper_w_p+lower_w_p);
                middle_w_f = 0.5*(upper_w_f+lower_w_f);
                initialise();
                rk4(middle_w_p,middle_w_f);
                fix();
                p_num_cross = pcrossings();
                f_num_cross = fcrossings();
                if (p_num_cross > 0)
                {
                       upper_w_p = middle_w_p;
                }
                else
                {
                       lower_w_p = middle_w_p;
                }
            // std::cout << "wp and wf : " << middle_w_p << " & "<< middle_w_f << std::endl;

            if (iter>65) break;
        }

        // std::cout << "#######################" << std::endl;
        // std::cout << "swapping over" << std::endl;
        // std::cout << "#######################" << std::endl;


        iter=0;
        lower_w_f = middle_w_f*0.5;
        upper_w_f = middle_w_f*2.;
        while (true)
        {
                iter++;
                middle_w_p = 0.5*(upper_w_p+lower_w_p);
                middle_w_f = 0.5*(upper_w_f+lower_w_f);
                initialise();
                rk4(middle_w_p,middle_w_f);
                fix();
                p_num_cross = pcrossings();
                f_num_cross = fcrossings();
                if (f_num_cross > 0)
                {
                       upper_w_f = middle_w_f;
                }
                else
                {
                       lower_w_f = middle_w_f;
                }
            // std::cout << "wp and wf : " << middle_w_p << " & "<< middle_w_f << std::endl;

            if (iter>65) break;
        }
    }


    WP = middle_w_p;
    WF = middle_w_f;
}


double QBSSolution::w_IB_soliton(double lower_w_, double upper_w_)
{
  int iter = 0, itermax, eigenstate, decimal_places_of_omega = 25.;
  double middle_w_;
  bool crossed;

  itermax = (int)((log(upper_w_)+decimal_places_of_omega*log(10.))/log(2.)); //calculate number of bisections needed (simple pen and paper calculation)
  while (true)
  {
      iter++;
      middle_w_ = 0.5*(upper_w_+lower_w_);
      initialise();
      rk4(middle_w_,-middle_w_);
      fix();
      crossed = soliton_eigen();
      // std::cout << "IB w : " << middle_w_ << ", " ;
      // std::cout << "crossed : " << crossed << std::endl;
      if (crossed)
      {
           upper_w_ = middle_w_;
      }
      else
      {
           lower_w_ = middle_w_;
      }
      if ((upper_w_-lower_w_)<w_tolerance) return upper_w_;
      if (upper_w_==lower_w_) return upper_w_;
      if (iter>100) return upper_w_;
  }
}

double QBSSolution::soliton_eigen()
{
    for (int i=1; i<gridsize; i++)
    {
        if (p[i]*p[i+1]<0.) return true;
        if (p[i]>p[i-1]) return false;
    }
}


// integrate the full ODE system from r=0 to dx*gridsize, this may (probably will) blow up at laarge raadius, but it is fixed by other functions aafterwards.
// it has smaller stepsize for the first (adaptive_buffer) steps.
void QBSSolution::rk4(const double w_p, const double w_f)
{
  	double k1, k2, k3, k4, q1, q2, q3, q4, w1, w2, w3, w4, x_=0., h = dx/2.;
    double e1, e2, e3, e4, f1, f2, f3, f4;
  	const double DX = dx;
    double DX_ = DX;
  	double o1, o2, o3, o4, s1, s2, s3, s4, r1, r2, r3, r4, t1, t2, t3, t4;
  	int index, jmax=0;
  	radius_array[0] = 0.;

  	for (int i = 1; i < gridsize; i++)
    {
    		
  		h = 0.5*dx;
        x_ = (i-1)*dx;


  		k1 = DX_*P_RHS(x_,p[i-1],f[i-1],z[i-1],dp[i-1],df[i-1],dz[i-1],psi[i-1],dpsi[i-1],omega[i-1],w_p,w_f);
        t1 = DX_*DP_RHS(x_,p[i-1],f[i-1],z[i-1],dp[i-1],df[i-1],dz[i-1],psi[i-1],dpsi[i-1],omega[i-1],w_p,w_f);
        e1 = DX_*F_RHS(x_,p[i-1],f[i-1],z[i-1],dp[i-1],df[i-1],dz[i-1],psi[i-1],dpsi[i-1],omega[i-1],w_p,w_f);
        f1 = DX_*DF_RHS(x_,p[i-1],f[i-1],z[i-1],dp[i-1],df[i-1],dz[i-1],psi[i-1],dpsi[i-1],omega[i-1],w_p,w_f);
  		q1 = DX_*Z_RHS(x_,p[i-1],f[i-1],z[i-1],dp[i-1],df[i-1],dz[i-1],psi[i-1],dpsi[i-1],omega[i-1],w_p,w_f);
        w1 = DX_*DZ_RHS(x_,p[i-1],f[i-1],z[i-1],dp[i-1],df[i-1],dz[i-1],psi[i-1],dpsi[i-1],omega[i-1],w_p,w_f);
  		o1 = DX_*OMEGA_RHS(x_,p[i-1],f[i-1],z[i-1],dp[i-1],df[i-1],dz[i-1],psi[i-1],dpsi[i-1],omega[i-1],w_p,w_f);
  		s1 = DX_*PSI_RHS(x_,p[i-1],f[i-1],z[i-1],dp[i-1],df[i-1],dz[i-1],psi[i-1],dpsi[i-1],omega[i-1],w_p,w_f);
  		r1 = DX_*DPSI_RHS(x_,p[i-1],f[i-1],z[i-1],dp[i-1],df[i-1],dz[i-1],psi[i-1],dpsi[i-1],omega[i-1],w_p,w_f);

  		k2 = DX_*P_RHS(x_ + h,p[i-1] + k1/2., f[i-1] + e1/2., z[i-1] + q1/2., dp[i-1] + t1/2., df[i-1] + f1/2., dz[i-1] + w1/2., psi[i-1] + s1/2.,dpsi[i-1] + r1/2.,omega[i-1] + o1/2.,w_p,w_f);
        t2 = DX_*DP_RHS(x_ + h,p[i-1] + k1/2., f[i-1] + e1/2., z[i-1] + q1/2., dp[i-1] + t1/2., df[i-1] + f1/2., dz[i-1] + w1/2., psi[i-1] + s1/2.,dpsi[i-1] + r1/2.,omega[i-1] + o1/2.,w_p,w_f);
        e2 = DX_*F_RHS(x_ + h,p[i-1] + k1/2., f[i-1] + e1/2., z[i-1] + q1/2., dp[i-1] + t1/2., df[i-1] + f1/2., dz[i-1] + w1/2., psi[i-1] + s1/2.,dpsi[i-1] + r1/2.,omega[i-1] + o1/2.,w_p,w_f);
        f2 = DX_*DF_RHS(x_ + h,p[i-1] + k1/2., f[i-1] + e1/2., z[i-1] + q1/2., dp[i-1] + t1/2., df[i-1] + f1/2., dz[i-1] + w1/2., psi[i-1] + s1/2.,dpsi[i-1] + r1/2.,omega[i-1] + o1/2.,w_p,w_f);
        q2 = DX_*Z_RHS(x_ + h,p[i-1] + k1/2., f[i-1] + e1/2., z[i-1] + q1/2., dp[i-1] + t1/2., df[i-1] + f1/2., dz[i-1] + w1/2., psi[i-1] + s1/2.,dpsi[i-1] + r1/2.,omega[i-1] + o1/2.,w_p,w_f);
  	    w2 = DX_*DZ_RHS(x_ + h,p[i-1] + k1/2., f[i-1] + e1/2., z[i-1] + q1/2., dp[i-1] + t1/2., df[i-1] + f1/2., dz[i-1] + w1/2., psi[i-1] + s1/2.,dpsi[i-1] + r1/2.,omega[i-1] + o1/2.,w_p,w_f);
        o2 = DX_*OMEGA_RHS(x_ + h,p[i-1] + k1/2., f[i-1] + e1/2., z[i-1] + q1/2., dp[i-1] + t1/2., df[i-1] + f1/2., dz[i-1] + w1/2., psi[i-1] + s1/2.,dpsi[i-1] + r1/2.,omega[i-1] + o1/2.,w_p,w_f);
  		s2 = DX_*PSI_RHS(x_ + h,p[i-1] + k1/2., f[i-1] + e1/2., z[i-1] + q1/2., dp[i-1] + t1/2., df[i-1] + f1/2., dz[i-1] + w1/2., psi[i-1] + s1/2.,dpsi[i-1] + r1/2.,omega[i-1] + o1/2.,w_p,w_f);
  		r2 = DX_*DPSI_RHS(x_ + h,p[i-1] + k1/2., f[i-1] + e1/2., z[i-1] + q1/2., dp[i-1] + t1/2., df[i-1] + f1/2., dz[i-1] + w1/2., psi[i-1] + s1/2.,dpsi[i-1] + r1/2.,omega[i-1] + o1/2.,w_p,w_f);

  		k3 = DX_*P_RHS(x_ + h,p[i-1] + k2/2., f[i-1] + e2/2., z[i-1] + q2/2., dp[i-1] + t2/2., df[i-1] + f2/2., dz[i-1] + w2/2., psi[i-1] + s2/2.,dpsi[i-1] + r2/2.,omega[i-1] + o2/2.,w_p,w_f);
        t3 = DX_*DP_RHS(x_ + h,p[i-1] + k2/2., f[i-1] + e2/2., z[i-1] + q2/2., dp[i-1] + t2/2., df[i-1] + f2/2., dz[i-1] + w2/2., psi[i-1] + s2/2.,dpsi[i-1] + r2/2.,omega[i-1] + o2/2.,w_p,w_f);
        e3 = DX_*F_RHS(x_ + h,p[i-1] + k2/2., f[i-1] + e2/2., z[i-1] + q2/2., dp[i-1] + t2/2., df[i-1] + f2/2., dz[i-1] + w2/2., psi[i-1] + s2/2.,dpsi[i-1] + r2/2.,omega[i-1] + o2/2.,w_p,w_f);
        f3 = DX_*DF_RHS(x_ + h,p[i-1] + k2/2., f[i-1] + e2/2., z[i-1] + q2/2., dp[i-1] + t2/2., df[i-1] + f2/2., dz[i-1] + w2/2., psi[i-1] + s2/2.,dpsi[i-1] + r2/2.,omega[i-1] + o2/2.,w_p,w_f);
  		q3 = DX_*Z_RHS(x_ + h,p[i-1] + k2/2., f[i-1] + e2/2., z[i-1] + q2/2., dp[i-1] + t2/2., df[i-1] + f2/2., dz[i-1] + w2/2., psi[i-1] + s2/2.,dpsi[i-1] + r2/2.,omega[i-1] + o2/2.,w_p,w_f);
        w3 = DX_*DZ_RHS(x_ + h,p[i-1] + k2/2., f[i-1] + e2/2., z[i-1] + q2/2., dp[i-1] + t2/2., df[i-1] + f2/2., dz[i-1] + w2/2., psi[i-1] + s2/2.,dpsi[i-1] + r2/2.,omega[i-1] + o2/2.,w_p,w_f);
  	    o3 = DX_*OMEGA_RHS(x_ + h,p[i-1] + k2/2., f[i-1] + e2/2., z[i-1] + q2/2., dp[i-1] + t2/2., df[i-1] + f2/2., dz[i-1] + w2/2., psi[i-1] + s2/2.,dpsi[i-1] + r2/2.,omega[i-1] + o2/2.,w_p,w_f);
  		s3 = DX_*PSI_RHS(x_ + h,p[i-1] + k2/2., f[i-1] + e2/2., z[i-1] + q2/2., dp[i-1] + t2/2., df[i-1] + f2/2., dz[i-1] + w2/2., psi[i-1] + s2/2.,dpsi[i-1] + r2/2.,omega[i-1] + o2/2.,w_p,w_f);
  		r3 = DX_*DPSI_RHS(x_ + h,p[i-1] + k2/2., f[i-1] + e2/2., z[i-1] + q2/2., dp[i-1] + t2/2., df[i-1] + f2/2., dz[i-1] + w2/2., psi[i-1] + s2/2.,dpsi[i-1] + r2/2.,omega[i-1] + o2/2.,w_p,w_f);

  		k4 = DX_*P_RHS(x_ + 2.*h,p[i-1] + k3, f[i-1] + e3, z[i-1] + q3, dp[i-1] + t3, df[i-1] + f3, dz[i-1] + w3, psi[i-1] + s3,dpsi[i-1] + r3,omega[i-1] + o3,w_p,w_f);
        t4 = DX_*DP_RHS(x_ + 2.*h,p[i-1] + k3, f[i-1] + e3, z[i-1] + q3, dp[i-1] + t3, df[i-1] + f3, dz[i-1] + w3, psi[i-1] + s3,dpsi[i-1] + r3,omega[i-1] + o3,w_p,w_f);
        e4 = DX_*F_RHS(x_ + 2.*h,p[i-1] + k3, f[i-1] + e3, z[i-1] + q3, dp[i-1] + t3, df[i-1] + f3, dz[i-1] + w3, psi[i-1] + s3,dpsi[i-1] + r3,omega[i-1] + o3,w_p,w_f);
        f4 = DX_*DF_RHS(x_ + 2.*h,p[i-1] + k3, f[i-1] + e3, z[i-1] + q3, dp[i-1] + t3, df[i-1] + f3, dz[i-1] + w3, psi[i-1] + s3,dpsi[i-1] + r3,omega[i-1] + o3,w_p,w_f);
  		q4 = DX_*Z_RHS(x_ + 2.*h,p[i-1] + k3, f[i-1] + e3, z[i-1] + q3, dp[i-1] + t3, df[i-1] + f3, dz[i-1] + w3, psi[i-1] + s3,dpsi[i-1] + r3,omega[i-1] + o3,w_p,w_f);
        w4 = DX_*DZ_RHS(x_ + 2.*h,p[i-1] + k3, f[i-1] + e3, z[i-1] + q3, dp[i-1] + t3, df[i-1] + f3, dz[i-1] + w3, psi[i-1] + s3,dpsi[i-1] + r3,omega[i-1] + o3,w_p,w_f);
  	    o4 = DX_*OMEGA_RHS(x_ + 2.*h,p[i-1] + k3, f[i-1] + e3, z[i-1] + q3, dp[i-1] + t3, df[i-1] + f3, dz[i-1] + w3, psi[i-1] + s3,dpsi[i-1] + r3,omega[i-1] + o3,w_p,w_f);
  		s4 = DX_*PSI_RHS(x_ + 2.*h,p[i-1] + k3, f[i-1] + e3, z[i-1] + q3, dp[i-1] + t3, df[i-1] + f3, dz[i-1] + w3, psi[i-1] + s3,dpsi[i-1] + r3,omega[i-1] + o3,w_p,w_f);
  		r4 = DX_*DPSI_RHS(x_ + 2.*h,p[i-1] + k3, f[i-1] + e3, z[i-1] + q3, dp[i-1] + t3, df[i-1] + f3, dz[i-1] + w3, psi[i-1] + s3,dpsi[i-1] + r3,omega[i-1] + o3,w_p,w_f);

  		p[i] = p[i-1] + (k1 + 2.*k2 + 2.*k3 + k4)/6.;
        dp[i] = dp[i-1] + (t1 + 2.*t2 + 2.*t3 + t4)/6.;
        f[i] = f[i-1] + (e1 + 2.*e2 + 2.*e3 + e4)/6.;
        df[i] = df[i-1] + (f1 + 2.*f2 + 2.*f3 + f4)/6.;
  		z[i] = z[i-1] + (q1 + 2.*q2 + 2.*q3 + q4)/6.;
        dz[i] = dz[i-1] + (w1 + 2.*w2 + 2.*w3 + w4)/6.;
  		psi[i] = psi[i-1] + (s1 + 2.*s2 + 2.*s3 + s4)/6.;
  		dpsi[i] = dpsi[i-1] + (r1 + 2.*r2 + 2.*r3 + r4)/6.;
  		omega[i] = omega[i-1] + (o1 + 2.*o2 + 2.*o3 + o4)/6.;
    
      	radius_array[i] = dx*i;
  	}
}


// same RK4 integrater, but scalar field is zero
void QBSSolution::rk4_VAC(const double w_p, const double w_f, const int i_start)
{
    double k1=0., k2=0., k3=0., k4=0., q1, q2, q3, q4, w1, w2, w3, w4, x_=0., h = dx/2.;
    const double DX = dx;
        double DX_ = DX;
    double o1, o2, o3, o4, s1, s2, s3, s4, r1, r2, r3, r4, t1=0., t2=0., t3=0., t4=0.;
    double e1=0., e2=0., e3=0., e4=0., f1=0., f2=0., f3=0., f4=0.;
    int index, jmax=0;
    radius_array[0] = 0.;

    //std::cout << "i_start : " << i_start << std::endl;
    for (int i = i_start; i < gridsize; i++)
    {
            
        h = 0.5*dx;
        x_ = (i-1)*dx;


        //k1 = DX_*P_RHS(x_,p[i-1],f[i-1],z[i-1],dp[i-1],df[i-1],dz[i-1],psi[i-1],dpsi[i-1],omega[i-1],w_p,w_f);
        //t1 = DX_*DP_RHS(x_,p[i-1],f[i-1],z[i-1],dp[i-1],df[i-1],dz[i-1],psi[i-1],dpsi[i-1],omega[i-1],w_p,w_f);
        //e1 = DX_*F_RHS(x_,p[i-1],f[i-1],z[i-1],dp[i-1],df[i-1],dz[i-1],psi[i-1],dpsi[i-1],omega[i-1],w_p,w_f);
        //f1 = DX_*DF_RHS(x_,p[i-1],f[i-1],z[i-1],dp[i-1],df[i-1],dz[i-1],psi[i-1],dpsi[i-1],omega[i-1],w_p,w_f);
        q1 = DX_*Z_RHS(x_,p[i-1],f[i-1],z[i-1],dp[i-1],df[i-1],dz[i-1],psi[i-1],dpsi[i-1],omega[i-1],w_p,w_f);
        w1 = DX_*DZ_RHS(x_,p[i-1],f[i-1],z[i-1],dp[i-1],df[i-1],dz[i-1],psi[i-1],dpsi[i-1],omega[i-1],w_p,w_f);
        o1 = DX_*OMEGA_RHS(x_,p[i-1],f[i-1],z[i-1],dp[i-1],df[i-1],dz[i-1],psi[i-1],dpsi[i-1],omega[i-1],w_p,w_f);
        s1 = DX_*PSI_RHS(x_,p[i-1],f[i-1],z[i-1],dp[i-1],df[i-1],dz[i-1],psi[i-1],dpsi[i-1],omega[i-1],w_p,w_f);
        r1 = DX_*DPSI_RHS(x_,p[i-1],f[i-1],z[i-1],dp[i-1],df[i-1],dz[i-1],psi[i-1],dpsi[i-1],omega[i-1],w_p,w_f);

        //k2 = DX_*P_RHS(x_ + h,p[i-1] + k1/2., f[i-1] + e1/2., z[i-1] + q1/2., dp[i-1] + t1/2., df[i-1] + f1/2., dz[i-1] + w1/2., psi[i-1] + s1/2.,dpsi[i-1] + r1/2.,omega[i-1] + o1/2.,w_p,w_f);
        //t2 = DX_*DP_RHS(x_ + h,p[i-1] + k1/2., f[i-1] + e1/2., z[i-1] + q1/2., dp[i-1] + t1/2., df[i-1] + f1/2., dz[i-1] + w1/2., psi[i-1] + s1/2.,dpsi[i-1] + r1/2.,omega[i-1] + o1/2.,w_p,w_f);
        //e2 = DX_*F_RHS(x_ + h,p[i-1] + k1/2., f[i-1] + e1/2., z[i-1] + q1/2., dp[i-1] + t1/2., df[i-1] + f1/2., dz[i-1] + w1/2., psi[i-1] + s1/2.,dpsi[i-1] + r1/2.,omega[i-1] + o1/2.,w_p,w_f);
        //f2 = DX_*DF_RHS(x_ + h,p[i-1] + k1/2., f[i-1] + e1/2., z[i-1] + q1/2., dp[i-1] + t1/2., df[i-1] + f1/2., dz[i-1] + w1/2., psi[i-1] + s1/2.,dpsi[i-1] + r1/2.,omega[i-1] + o1/2.,w_p,w_f);
        q2 = DX_*Z_RHS(x_ + h,p[i-1] + k1/2., f[i-1] + e1/2., z[i-1] + q1/2., dp[i-1] + t1/2., df[i-1] + f1/2., dz[i-1] + w1/2., psi[i-1] + s1/2.,dpsi[i-1] + r1/2.,omega[i-1] + o1/2.,w_p,w_f);
        w2 = DX_*DZ_RHS(x_ + h,p[i-1] + k1/2., f[i-1] + e1/2., z[i-1] + q1/2., dp[i-1] + t1/2., df[i-1] + f1/2., dz[i-1] + w1/2., psi[i-1] + s1/2.,dpsi[i-1] + r1/2.,omega[i-1] + o1/2.,w_p,w_f);
        o2 = DX_*OMEGA_RHS(x_ + h,p[i-1] + k1/2., f[i-1] + e1/2., z[i-1] + q1/2., dp[i-1] + t1/2., df[i-1] + f1/2., dz[i-1] + w1/2., psi[i-1] + s1/2.,dpsi[i-1] + r1/2.,omega[i-1] + o1/2.,w_p,w_f);
        s2 = DX_*PSI_RHS(x_ + h,p[i-1] + k1/2., f[i-1] + e1/2., z[i-1] + q1/2., dp[i-1] + t1/2., df[i-1] + f1/2., dz[i-1] + w1/2., psi[i-1] + s1/2.,dpsi[i-1] + r1/2.,omega[i-1] + o1/2.,w_p,w_f);
        r2 = DX_*DPSI_RHS(x_ + h,p[i-1] + k1/2., f[i-1] + e1/2., z[i-1] + q1/2., dp[i-1] + t1/2., df[i-1] + f1/2., dz[i-1] + w1/2., psi[i-1] + s1/2.,dpsi[i-1] + r1/2.,omega[i-1] + o1/2.,w_p,w_f);

        //k3 = DX_*P_RHS(x_ + h,p[i-1] + k2/2., f[i-1] + e2/2., z[i-1] + q2/2., dp[i-1] + t2/2., df[i-1] + f2/2., dz[i-1] + w2/2., psi[i-1] + s2/2.,dpsi[i-1] + r2/2.,omega[i-1] + o2/2.,w_p,w_f);
        //t3 = DX_*DP_RHS(x_ + h,p[i-1] + k2/2., f[i-1] + e2/2., z[i-1] + q2/2., dp[i-1] + t2/2., df[i-1] + f2/2., dz[i-1] + w2/2., psi[i-1] + s2/2.,dpsi[i-1] + r2/2.,omega[i-1] + o2/2.,w_p,w_f);
        //e3 = DX_*F_RHS(x_ + h,p[i-1] + k2/2., f[i-1] + e2/2., z[i-1] + q2/2., dp[i-1] + t2/2., df[i-1] + f2/2., dz[i-1] + w2/2., psi[i-1] + s2/2.,dpsi[i-1] + r2/2.,omega[i-1] + o2/2.,w_p,w_f);
        //f3 = DX_*DF_RHS(x_ + h,p[i-1] + k2/2., f[i-1] + e2/2., z[i-1] + q2/2., dp[i-1] + t2/2., df[i-1] + f2/2., dz[i-1] + w2/2., psi[i-1] + s2/2.,dpsi[i-1] + r2/2.,omega[i-1] + o2/2.,w_p,w_f);
        q3 = DX_*Z_RHS(x_ + h,p[i-1] + k2/2., f[i-1] + e2/2., z[i-1] + q2/2., dp[i-1] + t2/2., df[i-1] + f2/2., dz[i-1] + w2/2., psi[i-1] + s2/2.,dpsi[i-1] + r2/2.,omega[i-1] + o2/2.,w_p,w_f);
        w3 = DX_*DZ_RHS(x_ + h,p[i-1] + k2/2., f[i-1] + e2/2., z[i-1] + q2/2., dp[i-1] + t2/2., df[i-1] + f2/2., dz[i-1] + w2/2., psi[i-1] + s2/2.,dpsi[i-1] + r2/2.,omega[i-1] + o2/2.,w_p,w_f);
        o3 = DX_*OMEGA_RHS(x_ + h,p[i-1] + k2/2., f[i-1] + e2/2., z[i-1] + q2/2., dp[i-1] + t2/2., df[i-1] + f2/2., dz[i-1] + w2/2., psi[i-1] + s2/2.,dpsi[i-1] + r2/2.,omega[i-1] + o2/2.,w_p,w_f);
        s3 = DX_*PSI_RHS(x_ + h,p[i-1] + k2/2., f[i-1] + e2/2., z[i-1] + q2/2., dp[i-1] + t2/2., df[i-1] + f2/2., dz[i-1] + w2/2., psi[i-1] + s2/2.,dpsi[i-1] + r2/2.,omega[i-1] + o2/2.,w_p,w_f);
        r3 = DX_*DPSI_RHS(x_ + h,p[i-1] + k2/2., f[i-1] + e2/2., z[i-1] + q2/2., dp[i-1] + t2/2., df[i-1] + f2/2., dz[i-1] + w2/2., psi[i-1] + s2/2.,dpsi[i-1] + r2/2.,omega[i-1] + o2/2.,w_p,w_f);

        //k4 = DX_*P_RHS(x_ + 2.*h,p[i-1] + k3, f[i-1] + e3, z[i-1] + q3, dp[i-1] + t3, df[i-1] + f3, dz[i-1] + w3, psi[i-1] + s3,dpsi[i-1] + r3,omega[i-1] + o3,w_p,w_f);
        //t4 = DX_*DP_RHS(x_ + 2.*h,p[i-1] + k3, f[i-1] + e3, z[i-1] + q3, dp[i-1] + t3, df[i-1] + f3, dz[i-1] + w3, psi[i-1] + s3,dpsi[i-1] + r3,omega[i-1] + o3,w_p,w_f);
        //e4 = DX_*F_RHS(x_ + 2.*h,p[i-1] + k3, f[i-1] + e3, z[i-1] + q3, dp[i-1] + t3, df[i-1] + f3, dz[i-1] + w3, psi[i-1] + s3,dpsi[i-1] + r3,omega[i-1] + o3,w_p,w_f);
        //f4 = DX_*DF_RHS(x_ + 2.*h,p[i-1] + k3, f[i-1] + e3, z[i-1] + q3, dp[i-1] + t3, df[i-1] + f3, dz[i-1] + w3, psi[i-1] + s3,dpsi[i-1] + r3,omega[i-1] + o3,w_p,w_f);
        q4 = DX_*Z_RHS(x_ + 2.*h,p[i-1] + k3, f[i-1] + e3, z[i-1] + q3, dp[i-1] + t3, df[i-1] + f3, dz[i-1] + w3, psi[i-1] + s3,dpsi[i-1] + r3,omega[i-1] + o3,w_p,w_f);
        w4 = DX_*DZ_RHS(x_ + 2.*h,p[i-1] + k3, f[i-1] + e3, z[i-1] + q3, dp[i-1] + t3, df[i-1] + f3, dz[i-1] + w3, psi[i-1] + s3,dpsi[i-1] + r3,omega[i-1] + o3,w_p,w_f);
        o4 = DX_*OMEGA_RHS(x_ + 2.*h,p[i-1] + k3, f[i-1] + e3, z[i-1] + q3, dp[i-1] + t3, df[i-1] + f3, dz[i-1] + w3, psi[i-1] + s3,dpsi[i-1] + r3,omega[i-1] + o3,w_p,w_f);
        s4 = DX_*PSI_RHS(x_ + 2.*h,p[i-1] + k3, f[i-1] + e3, z[i-1] + q3, dp[i-1] + t3, df[i-1] + f3, dz[i-1] + w3, psi[i-1] + s3,dpsi[i-1] + r3,omega[i-1] + o3,w_p,w_f);
        r4 = DX_*DPSI_RHS(x_ + 2.*h,p[i-1] + k3, f[i-1] + e3, z[i-1] + q3, dp[i-1] + t3, df[i-1] + f3, dz[i-1] + w3, psi[i-1] + s3,dpsi[i-1] + r3,omega[i-1] + o3,w_p,w_f);

    
        radius_array[i] = dx*i;

        p[i] = 0.;//p[i-1] + (k1 + 2.*k2 + 2.*k3 + k4)/6.;
        dp[i] = 0.;//dp[i-1] + (t1 + 2.*t2 + 2.*t3 + t4)/6.;
        f[i] = 0.;//f[i-1] + (e1 + 2.*e2 + 2.*e3 + e4)/6.;
        df[i] = 0.;//df[i-1] + (f1 + 2.*f2 + 2.*f3 + f4)/6.;
        z[i] = z[i-1] + (q1 + 2.*q2 + 2.*q3 + q4)/6.;
        dz[i] = dz[i-1] + (w1 + 2.*w2 + 2.*w3 + w4)/6.;
        psi[i] = psi[i-1] + (s1 + 2.*s2 + 2.*s3 + s4)/6.;
        dpsi[i] = dpsi[i-1] + (r1 + 2.*r2 + 2.*r3 + r4)/6.;
        omega[i] = omega[i-1] + (o1 + 2.*o2 + 2.*o3 + o4)/6.;
    }

}



// takes an integrated ODE system and starts at point (iter) and re-integrates but enforcing vector field to decaay or be in vacuum
// the integral is adaptive in that it aaccelerates ar later radius in order to find correct asymptotic behaviour. It will shout if the radius reached is below 10e10
// bool adaaptive is true if stepsize is supposed to be adaptive aat large radius and false for constant stepsize
void QBSSolution::rk4_asymp(const double w_p, const double w_f, const int iter)
{
    double k1=0., k2=0., k3=0., k4=0., q1=0., q2=0., q3=0., q4=0., x_=iter*dx, h, delta = (double)gridsize;
    double t1=0., t2=0., t3=0., t4=0.;
    double e1=0., e2=0., e3=0., e4=0., f1=0., f2=0., f3=0., f4=0.;
    const double DX = dx;
    double DX_= DX;
    double o1, o2, o3, o4, s1, s2, s3, s4, r1, r2, r3, r4, w1, w2, w3, w4;
    double N_ = gridsize-iter, L_ = pow(9.,9);
    int i_;

    //std::cout << "In asymptotic calculator with start int " << iter << std::endl;

    double k_ = log(L_)/N_;

  	for (int i = iter+1; i < gridsize; ++i)
  	{
        i_ = double(i-iter);
	
		if (x_<9e9)
		{
			    DX_ = (exp(k_)-1.)*exp(k_*i_);
		}
		else
		{
			    DX_ = DX;
		}
		
		h = DX_/2.;


        //k1 = DX_*P_RHS(x_,p[i-1],f[i-1],z[i-1],dp[i-1],df[i-1],dz[i-1],psi[i-1],dpsi[i-1],omega[i-1],w_p,w_f);
        //t1 = DX_*DP_RHS(x_,p[i-1],f[i-1],z[i-1],dp[i-1],df[i-1],dz[i-1],psi[i-1],dpsi[i-1],omega[i-1],w_p,w_f);
        //e1 = DX_*F_RHS(x_,p[i-1],f[i-1],z[i-1],dp[i-1],df[i-1],dz[i-1],psi[i-1],dpsi[i-1],omega[i-1],w_p,w_f);
        //f1 = DX_*DF_RHS(x_,p[i-1],f[i-1],z[i-1],dp[i-1],df[i-1],dz[i-1],psi[i-1],dpsi[i-1],omega[i-1],w_p,w_f);
        q1 = DX_*Z_RHS(x_,p[i-1],f[i-1],z[i-1],dp[i-1],df[i-1],dz[i-1],psi[i-1],dpsi[i-1],omega[i-1],w_p,w_f);
        w1 = DX_*DZ_RHS(x_,p[i-1],f[i-1],z[i-1],dp[i-1],df[i-1],dz[i-1],psi[i-1],dpsi[i-1],omega[i-1],w_p,w_f);
        o1 = DX_*OMEGA_RHS(x_,p[i-1],f[i-1],z[i-1],dp[i-1],df[i-1],dz[i-1],psi[i-1],dpsi[i-1],omega[i-1],w_p,w_f);
        s1 = DX_*PSI_RHS(x_,p[i-1],f[i-1],z[i-1],dp[i-1],df[i-1],dz[i-1],psi[i-1],dpsi[i-1],omega[i-1],w_p,w_f);
        r1 = DX_*DPSI_RHS(x_,p[i-1],f[i-1],z[i-1],dp[i-1],df[i-1],dz[i-1],psi[i-1],dpsi[i-1],omega[i-1],w_p,w_f);

        //k2 = DX_*P_RHS(x_ + h,p[i-1] + k1/2., f[i-1] + e1/2., z[i-1] + q1/2., dp[i-1] + t1/2., df[i-1] + f1/2., dz[i-1] + w1/2., psi[i-1] + s1/2.,dpsi[i-1] + r1/2.,omega[i-1] + o1/2.,w_p,w_f);
        //t2 = DX_*DP_RHS(x_ + h,p[i-1] + k1/2., f[i-1] + e1/2., z[i-1] + q1/2., dp[i-1] + t1/2., df[i-1] + f1/2., dz[i-1] + w1/2., psi[i-1] + s1/2.,dpsi[i-1] + r1/2.,omega[i-1] + o1/2.,w_p,w_f);
        //e2 = DX_*F_RHS(x_ + h,p[i-1] + k1/2., f[i-1] + e1/2., z[i-1] + q1/2., dp[i-1] + t1/2., df[i-1] + f1/2., dz[i-1] + w1/2., psi[i-1] + s1/2.,dpsi[i-1] + r1/2.,omega[i-1] + o1/2.,w_p,w_f);
        //f2 = DX_*DF_RHS(x_ + h,p[i-1] + k1/2., f[i-1] + e1/2., z[i-1] + q1/2., dp[i-1] + t1/2., df[i-1] + f1/2., dz[i-1] + w1/2., psi[i-1] + s1/2.,dpsi[i-1] + r1/2.,omega[i-1] + o1/2.,w_p,w_f);
        q2 = DX_*Z_RHS(x_ + h,p[i-1] + k1/2., f[i-1] + e1/2., z[i-1] + q1/2., dp[i-1] + t1/2., df[i-1] + f1/2., dz[i-1] + w1/2., psi[i-1] + s1/2.,dpsi[i-1] + r1/2.,omega[i-1] + o1/2.,w_p,w_f);
        w2 = DX_*DZ_RHS(x_ + h,p[i-1] + k1/2., f[i-1] + e1/2., z[i-1] + q1/2., dp[i-1] + t1/2., df[i-1] + f1/2., dz[i-1] + w1/2., psi[i-1] + s1/2.,dpsi[i-1] + r1/2.,omega[i-1] + o1/2.,w_p,w_f);
        o2 = DX_*OMEGA_RHS(x_ + h,p[i-1] + k1/2., f[i-1] + e1/2., z[i-1] + q1/2., dp[i-1] + t1/2., df[i-1] + f1/2., dz[i-1] + w1/2., psi[i-1] + s1/2.,dpsi[i-1] + r1/2.,omega[i-1] + o1/2.,w_p,w_f);
        s2 = DX_*PSI_RHS(x_ + h,p[i-1] + k1/2., f[i-1] + e1/2., z[i-1] + q1/2., dp[i-1] + t1/2., df[i-1] + f1/2., dz[i-1] + w1/2., psi[i-1] + s1/2.,dpsi[i-1] + r1/2.,omega[i-1] + o1/2.,w_p,w_f);
        r2 = DX_*DPSI_RHS(x_ + h,p[i-1] + k1/2., f[i-1] + e1/2., z[i-1] + q1/2., dp[i-1] + t1/2., df[i-1] + f1/2., dz[i-1] + w1/2., psi[i-1] + s1/2.,dpsi[i-1] + r1/2.,omega[i-1] + o1/2.,w_p,w_f);

        //k3 = DX_*P_RHS(x_ + h,p[i-1] + k2/2., f[i-1] + e2/2., z[i-1] + q2/2., dp[i-1] + t2/2., df[i-1] + f2/2., dz[i-1] + w2/2., psi[i-1] + s2/2.,dpsi[i-1] + r2/2.,omega[i-1] + o2/2.,w_p,w_f);
        //t3 = DX_*DP_RHS(x_ + h,p[i-1] + k2/2., f[i-1] + e2/2., z[i-1] + q2/2., dp[i-1] + t2/2., df[i-1] + f2/2., dz[i-1] + w2/2., psi[i-1] + s2/2.,dpsi[i-1] + r2/2.,omega[i-1] + o2/2.,w_p,w_f);
        //e3 = DX_*F_RHS(x_ + h,p[i-1] + k2/2., f[i-1] + e2/2., z[i-1] + q2/2., dp[i-1] + t2/2., df[i-1] + f2/2., dz[i-1] + w2/2., psi[i-1] + s2/2.,dpsi[i-1] + r2/2.,omega[i-1] + o2/2.,w_p,w_f);
        //f3 = DX_*DF_RHS(x_ + h,p[i-1] + k2/2., f[i-1] + e2/2., z[i-1] + q2/2., dp[i-1] + t2/2., df[i-1] + f2/2., dz[i-1] + w2/2., psi[i-1] + s2/2.,dpsi[i-1] + r2/2.,omega[i-1] + o2/2.,w_p,w_f);
        q3 = DX_*Z_RHS(x_ + h,p[i-1] + k2/2., f[i-1] + e2/2., z[i-1] + q2/2., dp[i-1] + t2/2., df[i-1] + f2/2., dz[i-1] + w2/2., psi[i-1] + s2/2.,dpsi[i-1] + r2/2.,omega[i-1] + o2/2.,w_p,w_f);
        w3 = DX_*DZ_RHS(x_ + h,p[i-1] + k2/2., f[i-1] + e2/2., z[i-1] + q2/2., dp[i-1] + t2/2., df[i-1] + f2/2., dz[i-1] + w2/2., psi[i-1] + s2/2.,dpsi[i-1] + r2/2.,omega[i-1] + o2/2.,w_p,w_f);
        o3 = DX_*OMEGA_RHS(x_ + h,p[i-1] + k2/2., f[i-1] + e2/2., z[i-1] + q2/2., dp[i-1] + t2/2., df[i-1] + f2/2., dz[i-1] + w2/2., psi[i-1] + s2/2.,dpsi[i-1] + r2/2.,omega[i-1] + o2/2.,w_p,w_f);
        s3 = DX_*PSI_RHS(x_ + h,p[i-1] + k2/2., f[i-1] + e2/2., z[i-1] + q2/2., dp[i-1] + t2/2., df[i-1] + f2/2., dz[i-1] + w2/2., psi[i-1] + s2/2.,dpsi[i-1] + r2/2.,omega[i-1] + o2/2.,w_p,w_f);
        r3 = DX_*DPSI_RHS(x_ + h,p[i-1] + k2/2., f[i-1] + e2/2., z[i-1] + q2/2., dp[i-1] + t2/2., df[i-1] + f2/2., dz[i-1] + w2/2., psi[i-1] + s2/2.,dpsi[i-1] + r2/2.,omega[i-1] + o2/2.,w_p,w_f);

        //k4 = DX_*P_RHS(x_ + 2.*h,p[i-1] + k3, f[i-1] + e3, z[i-1] + q3, dp[i-1] + t3, df[i-1] + f3, dz[i-1] + w3, psi[i-1] + s3,dpsi[i-1] + r3,omega[i-1] + o3,w_p,w_f);
        //t4 = DX_*DP_RHS(x_ + 2.*h,p[i-1] + k3, f[i-1] + e3, z[i-1] + q3, dp[i-1] + t3, df[i-1] + f3, dz[i-1] + w3, psi[i-1] + s3,dpsi[i-1] + r3,omega[i-1] + o3,w_p,w_f);
        //e4 = DX_*F_RHS(x_ + 2.*h,p[i-1] + k3, f[i-1] + e3, z[i-1] + q3, dp[i-1] + t3, df[i-1] + f3, dz[i-1] + w3, psi[i-1] + s3,dpsi[i-1] + r3,omega[i-1] + o3,w_p,w_f);
        //f4 = DX_*DF_RHS(x_ + 2.*h,p[i-1] + k3, f[i-1] + e3, z[i-1] + q3, dp[i-1] + t3, df[i-1] + f3, dz[i-1] + w3, psi[i-1] + s3,dpsi[i-1] + r3,omega[i-1] + o3,w_p,w_f);
        q4 = DX_*Z_RHS(x_ + 2.*h,p[i-1] + k3, f[i-1] + e3, z[i-1] + q3, dp[i-1] + t3, df[i-1] + f3, dz[i-1] + w3, psi[i-1] + s3,dpsi[i-1] + r3,omega[i-1] + o3,w_p,w_f);
        w4 = DX_*DZ_RHS(x_ + 2.*h,p[i-1] + k3, f[i-1] + e3, z[i-1] + q3, dp[i-1] + t3, df[i-1] + f3, dz[i-1] + w3, psi[i-1] + s3,dpsi[i-1] + r3,omega[i-1] + o3,w_p,w_f);
        o4 = DX_*OMEGA_RHS(x_ + 2.*h,p[i-1] + k3, f[i-1] + e3, z[i-1] + q3, dp[i-1] + t3, df[i-1] + f3, dz[i-1] + w3, psi[i-1] + s3,dpsi[i-1] + r3,omega[i-1] + o3,w_p,w_f);
        s4 = DX_*PSI_RHS(x_ + 2.*h,p[i-1] + k3, f[i-1] + e3, z[i-1] + q3, dp[i-1] + t3, df[i-1] + f3, dz[i-1] + w3, psi[i-1] + s3,dpsi[i-1] + r3,omega[i-1] + o3,w_p,w_f);
        r4 = DX_*DPSI_RHS(x_ + 2.*h,p[i-1] + k3, f[i-1] + e3, z[i-1] + q3, dp[i-1] + t3, df[i-1] + f3, dz[i-1] + w3, psi[i-1] + s3,dpsi[i-1] + r3,omega[i-1] + o3,w_p,w_f);

    
        radius_array[i] = dx*i;

        p[i] = 0.;//p[i-1] + (k1 + 2.*k2 + 2.*k3 + k4)/6.;
        dp[i] = 0.;//dp[i-1] + (t1 + 2.*t2 + 2.*t3 + t4)/6.;
        f[i] = 0.;//f[i-1] + (e1 + 2.*e2 + 2.*e3 + e4)/6.;
        df[i] = 0.;//df[i-1] + (f1 + 2.*f2 + 2.*f3 + f4)/6.;
        z[i] = z[i-1] + (q1 + 2.*q2 + 2.*q3 + q4)/6.;
        dz[i] = dz[i-1] + (w1 + 2.*w2 + 2.*w3 + w4)/6.;
        psi[i] = psi[i-1] + (s1 + 2.*s2 + 2.*s3 + s4)/6.;
        dpsi[i] = dpsi[i-1] + (r1 + 2.*r2 + 2.*r3 + r4)/6.;
        omega[i] = omega[i-1] + (o1 + 2.*o2 + 2.*o3 + o4)/6.;

		x_ += DX_;

  	}

  	if (x_ < 8e7)
  	{
  	        std::cout << "\33[30;41m" << " Asymptotic Radius Too Small" << "\x1B[0m" << std::endl;
            std::cout << "Max radius " << x_ << " reached in asymptotic calculation! " << std::endl;
            std::cout << "Try increasing physical domain size or resolution! " << std::endl;
  	}
}



// these functions return the right hand side of the ode's

double QBSSolution::P_RHS(const double x, const double P, const double F, const double Z, const double DP, const double DF, const double DZ, const double PSI, const double DPSI, const double OM, const double w_p, const double w_f)
{
    return DP;
}

double QBSSolution::DP_RHS(const double x, const double P, const double F, const double Z, const double DP, const double DF, const double DZ, const double PSI, const double DPSI, const double OM, const double w_p, const double w_f)
{
    double r = ((x==0.)?eps:x);
    double DOM = OMEGA_RHS(x,P,F,Z,DP,DF,DZ,PSI,DPSI,OM,w_p,w_f);
    double RHS = -DP*(2./r + DOM/OM + DPSI/PSI);
    RHS += P*PSI*PSI*(DVDP(P,F) - pow(w_p + q_p*Z,2)/OM/OM);
    return RHS;
}

double QBSSolution::F_RHS(const double x, const double P, const double F, const double Z, const double DP, const double DF, const double DZ, const double PSI, const double DPSI, const double OM, const double w_p, const double w_f)
{
    return DF;
}

double QBSSolution::DF_RHS(const double x, const double P, const double F, const double Z, const double DP, const double DF, const double DZ, const double PSI, const double DPSI, const double OM, const double w_p, const double w_f)
{
    double r = ((x==0.)?eps:x);
    double DOM = OMEGA_RHS(x,P,F,Z,DP,DF,DZ,PSI,DPSI,OM,w_p,w_f);
    double RHS = -DF*(2./r + DOM/OM + DPSI/PSI);
    RHS += F*PSI*PSI*(DVDF(P,F) - pow(w_f + q_f*Z,2)/OM/OM);
    return RHS;
}

double QBSSolution::Z_RHS(const double x, const double P, const double F, const double Z, const double DP, const double DF, const double DZ, const double PSI, const double DPSI, const double OM, const double w_p, const double w_f)
{
    return DZ;
}

double QBSSolution::DZ_RHS(const double x, const double P, const double F, const double Z, const double DP, const double DF, const double DZ, const double PSI, const double DPSI, const double OM, const double w_p, const double w_f)
{
    double r = ((x==0.)?eps:x);
    double DOM = OMEGA_RHS(x,P,F,Z,DP,DF,DZ,PSI,DPSI,OM,w_p,w_f);
    double RHS = (DOM/OM - DPSI/PSI - 2./r)*DZ;
    RHS += 2.*P*P*PSI*PSI*q_p*(w_p+q_p*Z);
    RHS += 2.*F*F*PSI*PSI*q_f*(w_f+q_f*Z);
    return RHS;
}


double QBSSolution::PSI_RHS(const double x, const double P, const double F, const double Z, const double DP, const double DF, const double DZ, const double PSI, const double DPSI, const double OM, const double w_p, const double w_f)
{
  	double RHS = DPSI;
  	return RHS;
}

double QBSSolution::DPSI_RHS(const double x, const double P, const double F, const double Z, const double DP, const double DF, const double DZ, const double PSI, const double DPSI, const double OM, const double w_p, const double w_f)
{
    double r = ((x==0.)?eps:x);
    double Ttt =  P*P*pow(w_p+q_p*Z,2) + F*F*pow(w_f+q_f*Z,2);
    Ttt += OM*OM*(DP*DP + DF*DF)/(PSI*PSI);
    Ttt += OM*OM*V(P,F) + DZ*DZ/(2.*PSI*PSI);
    return -4.*M_PI*GNEWT*Ttt*PSI*PSI*PSI/(OM*OM) - 2.*DPSI/r + 0.5*DPSI*DPSI/PSI;
}


double QBSSolution::OMEGA_RHS(const double x, const double P, const double F, const double Z, const double DP, const double DF, const double DZ, const double PSI, const double DPSI, const double OM, const double w_p, const double w_f)
{
    double Trr =  ( P*P*pow(w_p+q_p*Z,2) + F*F*pow(w_f+q_f*Z,2) )*PSI*PSI/OM/OM;
    Trr += DP*DP + DF*DF - V(P,F)*PSI*PSI - DZ*DZ/(2.*OM*OM);
    return OM*(8.*M_PI*GNEWT*x*Trr*PSI*PSI - 2.*PSI*DPSI - x*DPSI*DPSI)/(2.*PSI*(x*DPSI+PSI));
}


// V is klein gordon potential and DV is its gradient. Depends on #define star_type at top
double QBSSolution::V(const double P, const double F)
{
    return m_p*m_p*P*P + m_f*m_f*F*F + lambda_fp*P*F*P*F + 0.5*(lambda_ff*F*F*F*F + lambda_pp*P*P*P*P);
}

// this is the derivativ dV/d(P^2)
double QBSSolution::DVDP(const double P, const double F)
{
    return m_p*m_p + lambda_fp*F*F + lambda_pp*P*P;
}

// this is the derivativ dV/d(F^2)
double QBSSolution::DVDF(const double P, const double F)
{
    return m_f*m_f + lambda_fp*P*P + lambda_ff*F*F;
}



// 4th order error (cubic interpolation) for field. shouts if asked to fetch a value outside the ode solution
double QBSSolution::get_p_interp(const double r) const
{
    int iter = (int) floor(r/dx); // index of 2nd (out of 4) gridpoints used for interpolation
    double a = (r/dx)-floor(r/dx)-0.5; //fraction from midpoint of two values, a = +- 1/2 is the nearest gridpoints
    double interpolated_value = 0, f1, f2, f3, f4;
    f1 = ((iter==0)?p[1]:p[iter-1]); // conditionl/ternary imposing zero gradeint at r=0
    f2 = p[iter];
    f3 = p[iter+1];
    f4 = p[iter+2];

    if (iter>gridsize-3){std::cout << "Trying to interpolate grid variable outside of initial data domain!" << std::endl;}

    // do the cubic spline, from mathematica script written by Robin (rc634@cam.ac.uk)
    interpolated_value = (1./48.)*(f1 *(-3.+2.*a+12.*a*a-8.*a*a*a) +(3.+2.*a)*(-(1.+2.*a)*(-9.*f3+f4+6.*f3*a-2*f4*a)+3.*f2*(3.-8.*a+4.*a*a)));
    return interpolated_value;
}

// 4th order error (cubic interpolation) for field. shouts if asked to fetch a value outside the ode solution
double QBSSolution::get_z_interp(const double r) const
{
    int iter = (int) floor(r/dx); // index of 2nd (out of 4) gridpoints used for interpolation
    double a = (r/dx)-floor(r/dx)-0.5; //fraction from midpoint of two values, a = +- 1/2 is the nearest gridpoints
    double interpolated_value = 0, f1, f2, f3, f4;
    f1 = ((iter==0)?z[1]:z[iter-1]); // conditionl/ternary imposing zero gradeint at r=0
    f2 = z[iter];
    f3 = z[iter+1];
    f4 = z[iter+2];

    if (iter>gridsize-3){std::cout << "Trying to interpolate grid variable outside of initial data domain!" << std::endl;}

    // do the cubic spline, from mathematica script written by Robin (rc634@cam.ac.uk)
    interpolated_value = (1./48.)*(f1 *(-3.+2.*a+12.*a*a-8.*a*a*a) +(3.+2.*a)*(-(1.+2.*a)*(-9.*f3+f4+6.*f3*a-2*f4*a)+3.*f2*(3.-8.*a+4.*a*a)));
    return interpolated_value;
}

double QBSSolution::get_lapse_interp(const double r) const
{
    int iter = (int) floor(r/dx); // index of 2nd (out of 4) gridpoints used for interpolation
    double a = (r/dx)-floor(r/dx)-0.5; //fraction from midpoint of two values, a = +- 1/2 is the nearest gridpoints
    double interpolated_value = 0, f1, f2, f3, f4;
    f1 = ((iter==0)?omega[1]:omega[iter-1]); // conditionl/ternary imposing zero gradeint at r=0
    f2 = omega[iter];
    f3 = omega[iter+1];
    f4 = omega[iter+2];

    if (iter>gridsize-3){std::cout << "Trying to interpolate grid variable outside of initial data domain!" << std::endl;}

    // do the cubic spline, from mathematica script written by Robin (rc634@cam.ac.uk)
    interpolated_value = (1./48.)*(f1 *(-3.+2.*a+12.*a*a-8.*a*a*a) +(3.+2.*a)*(-(1.+2.*a)*(-9.*f3+f4+6.*f3*a-2*f4*a)+3.*f2*(3.-8.*a+4.*a*a)));
    return interpolated_value;
}

double QBSSolution::get_psi_interp(const double r) const
{
    int iter = (int) floor(r/dx); // index of 2nd (out of 4) gridpoints used for interpolation
    double a = (r/dx)-floor(r/dx)-0.5; //fraction from midpoint of two values, a = +- 1/2 is the nearest gridpoints
    double interpolated_value = 0, f1, f2, f3, f4;
    f1 = ((iter==0)?psi[1]:psi[iter-1]); // conditionl/ternary imposing zero gradeint at r=0
    f2 = psi[iter];
    f3 = psi[iter+1];
    f4 = psi[iter+2];

    if (iter>gridsize-3){std::cout << "Trying to interpolate grid variable outside of initial data domain!" << std::endl;}

    // do the cubic spline, from mathematica script written by Robin (rc634@cam.ac.uk)
    interpolated_value = (1./48.)*(f1 *(-3.+2.*a+12.*a*a-8.*a*a*a) +(3.+2.*a)*(-(1.+2.*a)*(-9.*f3+f4+6.*f3*a-2*f4*a)+3.*f2*(3.-8.*a+4.*a*a)));
    return interpolated_value;
}




double QBSSolution::get_dpsi_interp(const double r) const
{
    int iter = (int) floor(r/dx); // index of 2nd (out of 4) gridpoints used for interpolation
    double a = (r/dx)-floor(r/dx)-0.5; //fraction from midpoint of two values, a = +- 1/2 is the nearest gridpoints
    double interpolated_value = 0, f1, f2, f3, f4;
    f1 = ((iter==0)?dpsi[1]:dpsi[iter-1]); // conditionl/ternary imposing zero gradeint at r=0
    f2 = dpsi[iter];
    f3 = dpsi[iter+1];
    f4 = dpsi[iter+2];

    if (iter>gridsize-3){std::cout << "Trying to interpolate grid variable outside of initial data domain!" << std::endl;}

    // do the cubic spline, from mathematica script written by Robin (rc634@cam.ac.uk)
    interpolated_value = (1./48.)*(f1 *(-3.+2.*a+12.*a*a-8.*a*a*a) +(3.+2.*a)*(-(1.+2.*a)*(-9.*f3+f4+6.*f3*a-2*f4*a)+3.*f2*(3.-8.*a+4.*a*a)));
    return interpolated_value;
}




double QBSSolution::get_dlapse_interp(const double r) const
{
    int iter = (int) floor(r/dx); // index of 2nd (out of 4) gridpoints used for interpolation
    double a = (r/dx)-floor(r/dx)-0.5; //fraction from midpoint of two values, a = +- 1/2 is the nearest gridpoints
    double interpolated_value = 0, f1, f2, f3, f4;
    f1 = ((iter==0)?omega[1]:omega[iter-1]);
    f2 = omega[iter];
    f3 = omega[iter+1];
    f4 = omega[iter+2];

    if (iter>gridsize-3){std::cout << "Trying to interpolate grid variable outside of initial data domain!" << std::endl;}

    // do the cubic spline (for gradient now), from mathematica script written by Robin (rc634@cam.ac.uk)
    interpolated_value = (1./(24.*dx))*( (f1-27.*f2+27.*f3-f4)  +  12.*a*(f1-f2-f3+f4)  -  12.*a*a*(f1-3.*f2+3.*f3-f4)  );
    return interpolated_value;
}



double QBSSolution::get_mass() const
{
    return ADM_MASS;
}

double QBSSolution::get_central_density() const
{
    return FINAL_FC;
}

double QBSSolution::get_w() const
{
    return w;
}

double QBSSolution::get_r(const double frac) const
{
    if ( (frac-0.5)*(frac-0.5) >= 0.25) {return -1.;}

    for (int i = 0; i < gridsize; ++i)
    {
        if (p[i]/p[0] < frac)
        {
            return radius_array[i];
        }
    }
}


void QBSSolution::set_initialcondition_params(QBSParams &m_params_QBS)
{
    gridsize = m_params_QBS.gridpoints;
    adaptive_buffer = 0.;//gridsize/10; // numer of gridpoints to intergate more carefully
    p.resize(gridsize); 
    dp.resize(gridsize); 
    f.resize(gridsize); 
    df.resize(gridsize); 
    z.resize(gridsize); 
    dz.resize(gridsize);
    psi.resize(gridsize); //conformal factor
    dpsi.resize(gridsize); //conformal factor gradient
    omega.resize(gridsize); // lapse
    radius_array.resize(gridsize); //radius
    radial_pressure.resize(gridsize); 
    N.resize(gridsize); 
    ham.resize(gridsize); 

    PC = m_params_QBS.p_central;
    FC = m_params_QBS.f_central;
    m_p = m_params_QBS.m_p;
    m_f = m_params_QBS.m_f;
    q_p = m_params_QBS.q_p;
    q_f = m_params_QBS.q_f;
    L = m_params_QBS.outer_radius; //just to make sure the function domain is slightly larger than the required cube
    dx = L/double((gridsize-1));
    lambda_ff = m_params_QBS.lambda_ff;
    lambda_fp = m_params_QBS.lambda_fp;
    lambda_pp = m_params_QBS.lambda_pp;
    repeats = m_params_QBS.solve_iterations;
}


void QBSSolution::shout() const
{
    std::cout << "Haliboombah!" << std::endl;
}


void QBSSolution::output_csv()
{
    std::ofstream p_file, q_file, f_file, psi_file, omega_file, r_file, press_file, ham_file, N_file;
    p_file.open("p.csv");
    f_file.open("f.csv");
    q_file.open("q.csv");
    psi_file.open("psi.csv");
    omega_file.open("omega.csv");
    r_file.open("r.csv");
    press_file.open("press.csv");
    ham_file.open("ham.csv");
    N_file.open("n.csv");

    for (int i=0; i<gridsize; i++)
    {
        p_file << p[i] << "," << std::endl;
        f_file << f[i] << "," << std::endl;
        q_file << z[i] << "," << std::endl;
        psi_file << psi[i] << "," << std::endl;
        omega_file << omega[i] << "," << std::endl;
        r_file << radius_array[i] << "," << std::endl;
        press_file << radial_pressure[i] << "," << std::endl;
        ham_file << ham[i] << "," << std::endl;
        N_file << N[i] << "," << std::endl;
    }

    p_file.close();
    f_file.close();
    q_file.close();
    psi_file.close();
    omega_file.close();
    r_file.close();
    press_file.close();
    ham_file.close();
    N_file.close();
}


#endif /* PROCASTARSOLUTION_IMPL_HPP_ */
