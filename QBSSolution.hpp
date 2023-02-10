
#ifndef QBSSOLUTION_HPP_
#define QBSSOLUTION_HPP_

#include <vector>
#include "cmath"


class QBSSolution
{

private: // private member variables/arrays
    double MM, FC, PC, ZC=0.0, q_p, q_f, m_p, m_f; //Klein Gordon mass squared, KG scalr field central amplitude
    const double GNEWT = 1;// 1./(M_PI*4.); // rescales KG field in ODE's
    double PSC=1., OMC=1.; // central density of scalar field (0.272 for kaup)  PSC and OMC are central values of conformal factor and lapse, not important as long as they are sensible (i.e. order 1)
    double lambda_ff, lambda_fp, lambda_pp; // phi 4 coupling in Klein gordon potential
    double sigma; // 0.2 works with PC = 0.05 // parameter for solitonic stars
    bool solitonic; // false fro mini/lambda star. true for solitonic star
    double EIGEN; // the desired eigenstate, 0 for ground
    int gridsize, adaptive_buffer, repeats; // anywhere from 2k-200k is ok
    const int adaptive_stepsize_repetitions = 20;//50; // 0 for no adaptive
  	double L, dx, W, w, WP=0., WF=0.; // L, length of domain, dx.
  	double OM_INF, PSI_INF, Z_INF; // asymptotics of lapse and cpnformal factpr
  	int mid_int; // integer where growing mode becomes relevant
  	double upper_w, lower_w, middle_w, upper_wp, lower_wp, middle_wp, upper_wf, lower_wf, middle_wf;
  	double adm_mass, aspect_mass;
    double eps = 10e-20;
    double w_tolerance = 10e-20;
    double fbork_minus_zbork;
    double ADM_MASS, FINAL_FC;


    std::vector<double> p; // field 
    std::vector<double> dp; // field deriv
    std::vector<double> f; // field 
    std::vector<double> df; // field deriv
    std::vector<double> z;  // leccy field
    std::vector<double> dz; // leccy field deriv
  	std::vector<double> N; // field 
    std::vector<double> mass_adm; // field 
    std::vector<double> ham; 
  	std::vector<double> psi; //conformal factor
  	std::vector<double> dpsi; //conformal factor gradient
  	std::vector<double> omega; // lapse
  	std::vector<double> radius_array; //radius


private: // private member fucntions functions
    void rk4(const double w_p, const double w_f);
    void rk4_VAC(const double w_p, const double w_f, const int i_start);
    void rk4_asymp(const double w_p, const double w_f, const int iter);
    

    double P_RHS(const double x, const double P, const double F, const double Z, const double DP, const double DF, const double DZ, const double PSI, const double DPSI, const double OM, const double w_p, const double w_f);
    double DP_RHS(const double x, const double P, const double F, const double Z, const double DP, const double DF, const double DZ, const double PSI, const double DPSI, const double OM, const double w_p, const double w_f);
    double F_RHS(const double x, const double P, const double F, const double Z, const double DP, const double DF, const double DZ, const double PSI, const double DPSI, const double OM, const double w_p, const double w_f);
    double DF_RHS(const double x, const double P, const double F, const double Z, const double DP, const double DF, const double DZ, const double PSI, const double DPSI, const double OM, const double w_p, const double w_f);
    double Z_RHS(const double x, const double P, const double F, const double Z, const double DP, const double DF, const double DZ, const double PSI, const double DPSI, const double OM, const double w_p, const double w_f);
    double DZ_RHS(const double x, const double P, const double F, const double Z, const double DP, const double DF, const double DZ, const double PSI, const double DPSI, const double OM, const double w_p, const double w_f);
    double PSI_RHS(const double x, const double P, const double F, const double Z, const double DP, const double DF, const double DZ, const double PSI, const double DPSI, const double OM, const double w_p, const double w_f);
    double DPSI_RHS(const double x, const double P, const double F, const double Z, const double DP, const double DF, const double DZ, const double PSI, const double DPSI, const double OM, const double w_p, const double w_f);
    double OMEGA_RHS(const double x, const double P, const double F, const double Z, const double DP, const double DF, const double DZ, const double PSI, const double DPSI, const double OM, const double w_p, const double w_f);


    void initialise();
    int pcrossings();
    int fcrossings();
    //bool f_crossed_g();
    void fix();
    void force_flat(const int iter_crit);
    void make_vacuum_outside_star(const double w_);
    double w_min(const double W_);
    double w_max(const double W_, const double lower_w_);
    double w_IB(double lower_w_, double upper_w_);
    double w_IB_bistar(const double repeats_param);
    double w_IB_bistar2();
    double w_IB_soliton(double lower_ww_, double upper_ww_);
    //double w_IB_f_equals_g(double lower_w_, double upper_w_,bool precise);
    double soliton_eigen();
    double find_W();
    double find_W_soliton();
    int find_midint();
    double V(const double P, const double F);
    double DVDP(const double P, const double F);
    double DVDF(const double P, const double F);


public:
    QBSSolution();
    void set_initialcondition_params(QBSParams &m_params_QBS);
    double get_p_interp(const double r) const;
    double get_z_interp(const double r) const;
    double get_lapse_interp(const double r) const;
    double get_psi_interp(const double r) const;
    double get_dpsi_interp(const double r) const;
    double get_dlapse_interp(const double r) const;
    double get_mass() const;
    double get_central_density() const;
    double get_w() const;
    double get_r(const double frac) const;
    void shout() const;
    void output_csv();
    void main();
};


#include "QBSSolution.impl.hpp"

#endif /* QBSSOLUTION_HPP_ */
