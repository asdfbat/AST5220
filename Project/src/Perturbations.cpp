#include <omp.h>
#include<string.h>
#include"Perturbations.h"

//====================================================
// Constructors
//====================================================

Perturbations::Perturbations(
    BackgroundCosmology *cosmo, 
    RecombinationHistory *rec) : 
  cosmo(cosmo), 
  rec(rec)
{}

//====================================================
// Do all the solving
//====================================================

void Perturbations::solve(){

  // Integrate all the perturbation equation and spline the result
  integrate_perturbations();

  // Compute source functions and spline the result
  // compute_source_functions();
}

//====================================================
// The main work: integrate all the perturbations
// and spline the results
//====================================================

void Perturbations::integrate_perturbations(){
  Utils::StartTiming("integrateperturbation");

  //===================================================================
  // TODO: Set up the k-array for the k's we are going to integrate over
  // Start at k_min end at k_max with n_k points with either a
  // quadratic or a logarithmic spacing
  //===================================================================
  Vector x_array_full = Utils::linspace(x_start, x_end, n_x);
  Vector log_k_array = Utils::linspace(log(k_min), log(k_max), n_k);
  Vector k_array(n_k);
  for(int i=0; i<n_k; i++)
    k_array.at(i) = exp(log_k_array.at(i));

  Vector2D delta_cdm_array(n_x, Vector(n_k));
  Vector2D delta_b_array(n_x, Vector(n_k));
  Vector2D v_cdm_array(n_x, Vector(n_k));
  Vector2D v_b_array(n_x, Vector(n_k));
  Vector2D Phi_array(n_x, Vector(n_k));
  Vector2D Psi_array(n_x, Vector(n_k));
  Vector2D Theta0_array(n_x, Vector(n_k));
  Vector2D Theta1_array(n_x, Vector(n_k));
  Vector2D Theta2_array(n_x, Vector(n_k));
  // Vector2D Theta3_array(n_x, Vector(n_k));
  // Vector2D Theta4_array(n_x, Vector(n_k));
  // Vector2D Theta5_array(n_x, Vector(n_k));
  // Vector2D Theta6_array(n_x, Vector(n_k));
  // Vector2D Theta7_array(n_x, Vector(n_k));
  // Vector3D Theta_array(Constants.n_ell_theta, Vector2D(n_x, Vector(n_k)));

  Vector delta_cdm_array_flat(n_x*n_k);
  Vector delta_b_array_flat(n_x*n_k);
  Vector v_cdm_array_flat(n_x*n_k);
  Vector v_b_array_flat(n_x*n_k);
  Vector Phi_array_flat(n_x*n_k);
  Vector Psi_array_flat(n_x*n_k);
  Vector Theta0_array_flat(n_x*n_k);
  Vector Theta1_array_flat(n_x*n_k);
  Vector Theta2_array_flat(n_x*n_k);
  // Vector Theta3_array_flat(n_x*n_k);
  // Vector Theta4_array_flat(n_x*n_k);
  // Vector Theta5_array_flat(n_x*n_k);
  // Vector Theta6_array_flat(n_x*n_k);
  // Vector Theta7_array_flat(n_x*n_k);
  // Vector2D Theta_array_flat(Constants.n_ell_theta, Vector(n_x*n_k));


  // Loop over all wavenumbers
  #pragma omp parallel for schedule(dynamic, 1)
  for(int ik = 0; ik < n_k; ik++){
    // Progress bar...
    if( (10*ik) / n_k != (10*ik+10) / n_k ) {
      std::cout << (100*ik+100)/n_k << "% " << std::flush;
      if(ik == n_k-1) std::cout << std::endl;
    }

    // Current value of k
    double k = k_array.at(ik);

    // Find value to integrate to
    double x_end_tight = get_tight_coupling_time(k);
    int idx_end_tight = (int) ((x_end_tight - x_start)/(x_end - x_start) * (n_x - 1.0));
    // printf("End of tight regime at x=%f, index=%d\n", x_end_tight, idx_end_tight);

    //===================================================================
    // TODO: Tight coupling integration
    // Remember to implement the routines:
    // set_ic : The IC at the start
    // rhs_tight_coupling_ode : The dydx for our coupled ODE system
    //===================================================================

    // Set up initial conditions in the tight coupling regime
    Vector y_tight_coupling_ini = set_ic(x_start, k);


    // The tight coupling ODE system
    ODEFunction dydx_tight_coupling = [&](double x, const double *y, double *dydx){
      return rhs_tight_coupling_ode(x, k, y, dydx);
    };

    // Integrate from x_start -> x_end_tight
    // ...
    // ...
    // ...
    // ...
    // ...
    Vector x_array_tc = Utils::linspace(x_start, x_end_tight, idx_end_tight+1);

    ODESolver ode_tc;
    ode_tc.solve(dydx_tight_coupling, x_array_tc, y_tight_coupling_ini);
    Vector tc_delta_cdm = ode_tc.get_data_by_component(Constants.ind_deltacdm_tc);
    Vector tc_delta_b   = ode_tc.get_data_by_component(Constants.ind_deltab_tc);
    Vector tc_v_cdm     = ode_tc.get_data_by_component(Constants.ind_vcdm_tc);
    Vector tc_v_b       = ode_tc.get_data_by_component(Constants.ind_vb_tc);
    Vector tc_Phi       = ode_tc.get_data_by_component(Constants.ind_Phi_tc);
    Vector tc_Theta0    = ode_tc.get_data_by_component(Constants.ind_start_theta_tc);
    Vector tc_Theta1    = ode_tc.get_data_by_component(Constants.ind_start_theta_tc + 1);

    double ck = Constants.c*k;

    for(int i=0; i<idx_end_tight+1; i++){
      double x = x_array_tc.at(i);
      double Hp = cosmo->Hp_of_x(x);
      double dtaudx = rec->dtaudx_of_x(x);

      delta_cdm_array.at(i).at(ik) = tc_delta_cdm.at(i);
      delta_b_array.at(i).at(ik)   = tc_delta_b.at(i);
      v_cdm_array.at(i).at(ik)     = tc_v_cdm.at(i);
      v_b_array.at(i).at(ik)       = tc_v_b.at(i);
      Phi_array.at(i).at(ik)       = tc_Phi.at(i);
      Psi_array.at(i).at(ik)       = -tc_Phi.at(i);
      Theta0_array.at(i).at(ik)    = tc_Theta0.at(i);
      Theta1_array.at(i).at(ik)    = tc_Theta1.at(i);
      Theta2_array.at(i).at(ik)    = -4.0*ck/(9.0*Hp*dtaudx)*tc_Theta1.at(i);
    }

    //===================================================================
    // TODO: Full equation integration
    // Remember to implement the routines:
    // set_ic_after_tight_coupling : The IC after tight coupling ends
    // rhs_full_ode : The dydx for our coupled ODE system
    //===================================================================

    // Set up initial conditions (y_tight_coupling is the solution at the end of tight coupling)
    // auto y_full_ini = set_ic_after_tight_coupling(y_tight_coupling, x_end_tight, k);

    // The full ODE system
    ODEFunction dydx_nontc = [&](double x, const double *y, double *dydx){
      return rhs_full_ode(x, k, y, dydx);
    };

    // Integrate from x_end_tight -> x_end
    // ...
    // ...
    // ...
    // ...
    // ...
    double x = x_array_full[idx_end_tight];
    double Hp = cosmo->Hp_of_x(x);
    double dtaudx = rec->dtaudx_of_x(x);

    // printf("%d\n", Constants.n_ell_tot_full);
    Vector y_nontc_ini(Constants.n_ell_tot_full);
    y_nontc_ini[Constants.ind_vcdm]          = tc_v_cdm[idx_end_tight];
    y_nontc_ini[Constants.ind_vb]            = tc_v_b[idx_end_tight];
    y_nontc_ini[Constants.ind_deltacdm]      = tc_delta_cdm[idx_end_tight];
    y_nontc_ini[Constants.ind_deltab]        = tc_delta_b[idx_end_tight];
    y_nontc_ini[Constants.ind_Phi]           = tc_Phi[idx_end_tight];
    y_nontc_ini[Constants.ind_start_theta]   = tc_Theta0[idx_end_tight];
    y_nontc_ini[Constants.ind_start_theta+1] = tc_Theta1[idx_end_tight];
    y_nontc_ini[Constants.ind_start_theta+2] = -4.0/9.0*ck/(Hp*dtaudx)*tc_Theta1[idx_end_tight];
    for(int l=3; l<Constants.n_ell_theta; l++){
      y_nontc_ini[Constants.ind_start_theta+l] = l/(2*l + 1)*ck/(Hp*dtaudx)*y_nontc_ini[Constants.ind_start_theta+l-1];
    }

    // printf("ASDF -> %lu\n", y_nontc_ini.size());
  

    int n_x_nontc = n_x - idx_end_tight;


    // printf("---- %d %d %d\n", n_x_nontc, n_x, idx_end_tight);


    Vector x_array_nontc = Utils::linspace(x_end_tight, x_end, n_x_nontc);
    ODESolver ode_nontc;
    ode_nontc.solve(dydx_nontc, x_array_nontc, y_nontc_ini);

    Vector nontc_delta_cdm = ode_nontc.get_data_by_component(Constants.ind_deltacdm);
    Vector nontc_delta_b   = ode_nontc.get_data_by_component(Constants.ind_deltab);
    Vector nontc_v_cdm     = ode_nontc.get_data_by_component(Constants.ind_vcdm);
    Vector nontc_v_b       = ode_nontc.get_data_by_component(Constants.ind_vb);
    Vector nontc_Phi       = ode_nontc.get_data_by_component(Constants.ind_Phi);
    Vector nontc_Theta0    = ode_nontc.get_data_by_component(Constants.ind_start_theta);
    Vector nontc_Theta1    = ode_nontc.get_data_by_component(Constants.ind_start_theta + 1);
    Vector nontc_Theta2    = ode_nontc.get_data_by_component(Constants.ind_start_theta + 2);
    // Vector nontc_Theta3    = ode_nontc.get_data_by_component(Constants.ind_start_theta + 3);
    // Vector nontc_Theta4    = ode_nontc.get_data_by_component(Constants.ind_start_theta + 4);
    // Vector nontc_Theta5    = ode_nontc.get_data_by_component(Constants.ind_start_theta + 5);
    // Vector nontc_Theta6    = ode_nontc.get_data_by_component(Constants.ind_start_theta + 6);
    // Vector nontc_Theta7    = ode_nontc.get_data_by_component(Constants.ind_start_theta + 7);


    for(int i=idx_end_tight+1; i<n_x; i++){
      // printf("\n%d ", i-idx_end_tight);
      delta_cdm_array.at(i).at(ik) = nontc_delta_cdm.at(i-idx_end_tight);
      delta_b_array.at(i).at(ik)   = nontc_delta_b.at(i-idx_end_tight);
      v_cdm_array.at(i).at(ik)     = nontc_v_cdm.at(i-idx_end_tight);
      v_b_array.at(i).at(ik)       = nontc_v_b.at(i-idx_end_tight);
      Phi_array.at(i).at(ik)       = nontc_Phi.at(i-idx_end_tight);
      Psi_array.at(i).at(ik)       = -nontc_Phi.at(i-idx_end_tight);
      Theta0_array.at(i).at(ik)    = nontc_Theta0.at(i-idx_end_tight);
      Theta1_array.at(i).at(ik)    = nontc_Theta1.at(i-idx_end_tight);
      Theta2_array.at(i).at(ik)    = nontc_Theta2.at(i-idx_end_tight);
      // Theta3_array.at(i).at(ik)    = nontc_Theta3.at(i-idx_end_tight-1);
      // Theta4_array.at(i).at(ik)    = nontc_Theta4.at(i-idx_end_tight-1);
      // Theta5_array.at(i).at(ik)    = nontc_Theta5.at(i-idx_end_tight-1);
      // Theta6_array.at(i).at(ik)    = nontc_Theta6.at(i-idx_end_tight-1);
      // Theta7_array.at(i).at(ik)    = nontc_Theta7.at(i-idx_end_tight-1);
    }


    //===================================================================
    // TODO: remember to store the data found from integrating so we can
    // spline it below
    //
    // To compute a 2D spline of a function f(x,k) the data must be given 
    // to the spline routine as a 1D array f_array with the points f(ix, ik) 
    // stored as f_array[ix + n_x * ik]
    // Example:
    // Vector x_array(n_x);
    // Vector k_array(n_k);
    // Vector f(n_x * n_k);
    // Spline2D y_spline;
    // f_spline.create(x_array, k_array, f_array);
    // We can now use the spline as f_spline(x, k)
    //
    //===================================================================
    //...
    //...


  }
  Utils::EndTiming("integrateperturbation");

  //=============================================================================
  // TODO: Make all splines needed: Theta0,Theta1,Theta2,Phi,Psi,...
  //=============================================================================
  // ...
  // ...
  // ...

  for(int j=0; j<n_k; j++){
    for(int i=0; i<n_x; i++){
      delta_cdm_array_flat.at(j*n_x + i) = delta_cdm_array.at(i).at(j);
      delta_b_array_flat.at(j*n_x + i) = delta_b_array.at(i).at(j);
      v_cdm_array_flat.at(j*n_x + i) = v_cdm_array.at(i).at(j);
      v_b_array_flat.at(j*n_x + i) = v_b_array.at(i).at(j);
      Phi_array_flat.at(j*n_x + i) = Phi_array.at(i).at(j);
      Psi_array_flat.at(j*n_x + i) = Psi_array.at(i).at(j);
      Theta0_array_flat.at(j*n_x + i) = Theta0_array.at(i).at(j);
      Theta1_array_flat.at(j*n_x + i) = Theta1_array.at(i).at(j);
      Theta2_array_flat.at(j*n_x + i) = Theta2_array.at(i).at(j);
      // Theta3_array_flat.at(j*n_x + i) = Theta3_array.at(i).at(j);
      // Theta4_array_flat.at(j*n_x + i) = Theta4_array.at(i).at(j);
      // Theta5_array_flat.at(j*n_x + i) = Theta5_array.at(i).at(j);
      // Theta6_array_flat.at(j*n_x + i) = Theta6_array.at(i).at(j);
      // Theta7_array_flat.at(j*n_x + i) = Theta7_array.at(i).at(j);
    }
  }


  delta_cdm_spline.create(x_array_full, log_k_array, delta_cdm_array_flat, "delta_cdm_spline");
  delta_b_spline.create(x_array_full, log_k_array, delta_b_array_flat, "delta_b_spline");
  v_cdm_spline.create(x_array_full, log_k_array, v_cdm_array_flat, "v_cdm_spline");
  v_b_spline.create(x_array_full, log_k_array, v_b_array_flat, "v_b_spline");
  Phi_spline.create(x_array_full, log_k_array, Phi_array_flat, "Phi_spline");
  Psi_spline.create(x_array_full, log_k_array, Psi_array_flat, "Psi_spline");
  Pi_spline.create(x_array_full, log_k_array, Theta2_array_flat, "Psi_spline");
  Theta0_spline.create(x_array_full, log_k_array, Theta0_array_flat, "Theta0_spline");
  Theta1_spline.create(x_array_full, log_k_array, Theta1_array_flat, "Theta1_spline");
  Theta2_spline.create(x_array_full, log_k_array, Theta2_array_flat, "Theta2_spline");
}

//====================================================
// Set IC at the start of the run (this is in the
// tight coupling regime)
//====================================================
Vector Perturbations::set_ic(const double x, const double k) const{

  // The vector we are going to fill
  Vector y_tc(Constants.n_ell_tot_tc);

  //=============================================================================
  // Compute where in the y_tc array each component belongs
  // This is just an example of how to do it to make it easier
  // Feel free to organize the component any way you like
  //=============================================================================
  // For integration of perturbations in tight coupling regime (Only 2 photon multipoles + neutrinos needed)
  const int n_ell_theta_tc      = Constants.n_ell_theta_tc;
  const int n_ell_neutrinos_tc  = Constants.n_ell_neutrinos_tc;
  const int n_ell_tot_tc        = Constants.n_ell_tot_tc;
  const bool polarization       = Constants.polarization;
  const bool neutrinos          = Constants.neutrinos;

  // References to the tight coupling quantities
  double &delta_cdm    =  y_tc[Constants.ind_deltacdm_tc];
  double &delta_b      =  y_tc[Constants.ind_deltab_tc];
  double &v_cdm        =  y_tc[Constants.ind_vcdm_tc];
  double &v_b          =  y_tc[Constants.ind_vb_tc];
  double &Phi          =  y_tc[Constants.ind_Phi_tc];
  double *Theta        = &y_tc[Constants.ind_start_theta_tc];

  //=============================================================================
  // TODO: Set the initial conditions in the tight coupling regime
  //=============================================================================
  // ...
  // ...
  const double ck     = Constants.c*k;
  const double dtaudx = rec->dtaudx_of_x(x);
  const double Hp     = cosmo->Hp_of_x(x);
  const double fv     = 0;
  const double Psi    = -1.0/(1.5 + 0.4*fv);

  Phi        = -(1.0 + 0.4*fv)*Psi;
  delta_b  = -1.5*Psi;
  delta_cdm  = delta_b;
  v_b      = -ck*Psi/(2.0*Hp);
  v_cdm      = v_b;
  Theta[0]   = -0.5*Psi;
  Theta[1]   = ck*Psi/(6.0*Hp);

  // printf("IC:\nk=%e\nTheta0=%e\nTheta1=%e\nDcdm=%e\nvcdm=%e\nPhi=%e\n", k*Constants.Mpc, Theta[0], Theta[1], delta_cdm, v_cdm, Phi);

  

  // SET: Scalar quantities (Gravitational potential, baryons and CDM)
  // ...
  // ...

  // SET: Photon temperature perturbations (Theta_ell)
  // ...
  // ...
  return y_tc;
}

//====================================================
// Set IC for the full ODE system after tight coupling 
// regime ends
//====================================================

Vector Perturbations::set_ic_after_tight_coupling(
    const Vector &y_tc, 
    const double x, 
    const double k) const{

  // Make the vector we are going to fill
  Vector y(Constants.n_ell_tot_full);
  
  //=============================================================================
  // Compute where in the y array each component belongs and where corresponding
  // components are located in the y_tc array
  // This is just an example of how to do it to make it easier
  // Feel free to organize the component any way you like
  //=============================================================================

  // Number of multipoles we have in the full regime
  const int n_ell_theta         = Constants.n_ell_theta;
  const int n_ell_thetap        = Constants.n_ell_thetap;
  const int n_ell_neutrinos     = Constants.n_ell_neutrinos;
  const bool polarization       = Constants.polarization;
  const bool neutrinos          = Constants.neutrinos;

  // Number of multipoles we have in the tight coupling regime
  const int n_ell_theta_tc      = Constants.n_ell_theta_tc;
  const int n_ell_neutrinos_tc  = Constants.n_ell_neutrinos_tc;

  // References to the tight coupling quantities
  const double &delta_cdm_tc    =  y_tc[Constants.ind_deltacdm_tc];
  const double &delta_b_tc      =  y_tc[Constants.ind_deltab_tc];
  const double &v_cdm_tc        =  y_tc[Constants.ind_vcdm_tc];
  const double &v_b_tc          =  y_tc[Constants.ind_vb_tc];
  const double &Phi_tc          =  y_tc[Constants.ind_Phi_tc];
  const double *Theta_tc        = &y_tc[Constants.ind_start_theta_tc];

  // References to the quantities we are going to set
  double &delta_cdm       =  y[Constants.ind_deltacdm_tc];
  double &delta_b         =  y[Constants.ind_deltab_tc];
  double &v_cdm           =  y[Constants.ind_vcdm_tc];
  double &v_b             =  y[Constants.ind_vb_tc];
  double &Phi             =  y[Constants.ind_Phi_tc];
  double *Theta           = &y[Constants.ind_start_theta_tc];

  //=============================================================================
  // TODO: fill in the initial conditions for the full equation system below
  // NB: remember that we have different number of multipoles in the two
  // regimes so be careful when assigning from the tc array
  //=============================================================================
  // ...
  // ...
  // ...

  // SET: Scalar quantities (Gravitational potental, baryons and CDM)
  // ...
  // ...

  // SET: Photon temperature perturbations (Theta_ell)
  // ...
  // ...

 

  return y;
}

//====================================================
// The time when tight coupling end
//====================================================

double Perturbations::get_tight_coupling_time(const double k) const{
  double x_tight_coupling_end = 0.0;

  //=============================================================================
  // TODO: compute and return x for when tight coupling ends
  // Remember all the three conditions in Callin
  //=============================================================================
  // ...
  // ...
  int npts = 1000;
  Vector x_arr = Utils::linspace(x_start, 0, npts);
  for(int i=0; i<npts; i++){
    double x = x_arr[i];
    double cond1 = abs(rec->dtaudx_of_x(x)/10.0);
    double cond2 = abs(rec->dtaudx_of_x(x)*cosmo->Hp_of_x(x)/(k*Constants.c*10.0));
    double cond3 = rec->Xe_of_x(x)/0.99;

    if((cond1 < 1) || (cond2 < 1) || (cond3 < 1)){
      return x;
    }
  }
  printf("ERROR: Could not find condition for end of tight coupling regime.");
  return -9999999;
  // return x_tight_coupling_end;
}

//====================================================
// After integrsating the perturbation compute the
// source function(s)
//====================================================
void Perturbations::compute_source_functions(){
  Utils::StartTiming("source");

  //=============================================================================
  // TODO: Make the x and k arrays to evaluate over and use to make the splines
  //=============================================================================
  // ...
  // ...
  Vector k_array;
  Vector x_array;

  // Make storage for the source functions (in 1D array to be able to pass it to the spline)
  Vector ST_array(k_array.size() * x_array.size());
  Vector SE_array(k_array.size() * x_array.size());

  // Compute source functions
  for(auto ix = 0; ix < x_array.size(); ix++){
    const double x = x_array[ix];
    for(auto ik = 0; ik < k_array.size(); ik++){
      const double k = k_array[ik];

      // NB: This is the format the data needs to be stored 
      // in a 1D array for the 2D spline routine source(ix,ik) -> S_array[ix + nx * ik]
      const int index = ix + n_x * ik;

      //=============================================================================
      // TODO: Compute the source functions
      //=============================================================================
      // Fetch all the things we need...
      // const double Hp       = cosmo->Hp_of_x(x);
      // const double tau      = rec->tau_of_x(x);
      // ...
      // ...

      // Temperatur source
      ST_array[index] = 0.0;

      // Polarization source
      if(Constants.polarization){
        SE_array[index] = 0.0;
      }
    }
  }

  // Spline the source functions
  ST_spline.create (x_array, k_array, ST_array, "Source_Temp_x_k");
  if(Constants.polarization){
    SE_spline.create (x_array, k_array, SE_array, "Source_Pol_x_k");
  }

  Utils::EndTiming("source");
}

//====================================================
// The right hand side of the perturbations ODE
// in the tight coupling regime
//====================================================

// Derivatives in the tight coupling regime
int Perturbations::rhs_tight_coupling_ode(double x, double k, const double *y, double *dydx){

  //=============================================================================
  // Compute where in the y / dydx array each component belongs
  // This is just an example of how to do it to make it easier
  // Feel free to organize the component any way you like
  //=============================================================================
  
  // For integration of perturbations in tight coupling regime (Only 2 photon multipoles + neutrinos needed)
  const int n_ell_theta_tc      = Constants.n_ell_theta_tc;
  const int n_ell_neutrinos_tc  = Constants.n_ell_neutrinos_tc;
  const bool neutrinos          = Constants.neutrinos;

  // The different quantities in the y array
  const double &delta_cdm       =  y[Constants.ind_deltacdm_tc];
  const double &delta_b         =  y[Constants.ind_deltab_tc];
  const double &v_cdm           =  y[Constants.ind_vcdm_tc];
  const double &v_b             =  y[Constants.ind_vb_tc];
  const double &Phi             =  y[Constants.ind_Phi_tc];
  const double *Theta           = &y[Constants.ind_start_theta_tc];

  // References to the quantities we are going to set in the dydx array
  double &ddelta_cdmdx    =  dydx[Constants.ind_deltacdm_tc];
  double &ddelta_bdx      =  dydx[Constants.ind_deltab_tc];
  double &dv_cdmdx        =  dydx[Constants.ind_vcdm_tc];
  double &dv_bdx          =  dydx[Constants.ind_vb_tc];
  double &dPhidx          =  dydx[Constants.ind_Phi_tc];
  double *dThetadx        = &dydx[Constants.ind_start_theta_tc];

  //=============================================================================
  // TODO: fill in the expressions for all the derivatives
  //=============================================================================
  const double ck = Constants.c*k;
  const double a  = exp(x);
  const double Hp = cosmo->Hp_of_x(x);
  const double dHpdx = cosmo->dHpdx_of_x(x);
  const double H0 = cosmo->get_H0();
  const double OmegaR   = cosmo->get_OmegaR();  //OMEGA0 OR OMEGA(x)???
  const double OmegaCDM = cosmo->get_OmegaCDM();
  const double OmegaB   = cosmo->get_OmegaB();
  const double dtaudx   = rec->dtaudx_of_x(x);
  const double ddtauddx = rec->ddtauddx_of_x(x);
  const double R  = 4*OmegaR/(3*OmegaB*a);


  // SET: Scalar quantities (Phi, delta, v, ...)
  // ...
  // ...
  // ...
  const double Theta2 = -4.0*ck/(9.0*Hp*dtaudx)*Theta[1];  // Theta[2] doesn't exist in tc array, so we set it here.
  const double Psi    = -Phi - 12.0*H0*H0/(ck*ck*a*a)*OmegaR*Theta2;

  dPhidx        = Psi - ck*ck/(3.0*Hp*Hp)*Phi + H0*H0/(2.0*Hp*Hp)*(OmegaCDM/a*delta_cdm + OmegaB/a*delta_b + 4.0*OmegaR/(a*a)*Theta[0]);
  ddelta_cdmdx  = ck*v_cdm/(Hp) - 3.0*dPhidx;
  dv_cdmdx      = -v_cdm - ck/Hp*Psi;
  ddelta_bdx    = ck/Hp*v_b - 3.0*dPhidx;
  dThetadx[0]   = -ck/Hp*Theta[1] - dPhidx;

  const double q = (-((1 - R)*dtaudx + (1 + R)*ddtauddx)*(3*Theta[1] + v_b) - ck/Hp*Psi + (1 - dHpdx/Hp)*ck/Hp*(-Theta[0] + 2*Theta2) - ck/Hp*dThetadx[0])/((1 + R)*dtaudx + dHpdx/Hp - 1);
  dv_bdx = 1/(1 + R)*(-v_b - ck/Hp*Psi + R*(q + ck/Hp*(-Theta[0] + 2*Theta2) - ck/Hp*Psi));
  dThetadx[1] = 1.0/3.0*(q - dv_bdx);


  return GSL_SUCCESS;
}

//====================================================
// The right hand side of the full ODE
//====================================================

int Perturbations::rhs_full_ode(double x, double k, const double *y, double *dydx){
  
  //=============================================================================
  // Compute where in the y / dydx array each component belongs
  // This is just an example of how to do it to make it easier
  // Feel free to organize the component any way you like
  //=============================================================================

  // Index and number of the different quantities
  const int n_ell_theta         = Constants.n_ell_theta;
  const int n_ell_thetap        = Constants.n_ell_thetap;
  const int n_ell_neutrinos     = Constants.n_ell_neutrinos;
  const bool polarization       = Constants.polarization;
  const bool neutrinos          = Constants.neutrinos;

  // The different quantities in the y array
  const double &delta_cdm       =  y[Constants.ind_deltacdm];
  const double &delta_b         =  y[Constants.ind_deltab];
  const double &v_cdm           =  y[Constants.ind_vcdm];
  const double &v_b             =  y[Constants.ind_vb];
  const double &Phi             =  y[Constants.ind_Phi];
  const double *Theta           = &y[Constants.ind_start_theta];

  // References to the quantities we are going to set in the dydx array
  double &ddelta_cdmdx    =  dydx[Constants.ind_deltacdm];
  double &ddelta_bdx      =  dydx[Constants.ind_deltab];
  double &dv_cdmdx        =  dydx[Constants.ind_vcdm];
  double &dv_bdx          =  dydx[Constants.ind_vb];
  double &dPhidx          =  dydx[Constants.ind_Phi];
  double *dThetadx        = &dydx[Constants.ind_start_theta];

  // Cosmological parameters and variables
  // double Hp = cosmo->Hp_of_x(x);
  // ...

  // Recombination variables
  // ...

  //=============================================================================
  // TODO: fill in the expressions for all the derivatives
  //=============================================================================

  // SET: Scalar quantities (Phi, delta, v, ...)
  // ...
  // ...
  // ...

  // SET: Photon multipoles (Theta_ell)
  // ...
  // ...
  // ...

  const double ck = Constants.c*k;
  const double a  = exp(x);
  const double Hp = cosmo->Hp_of_x(x);
  const double dHpdx = cosmo->dHpdx_of_x(x);
  const double H0 = cosmo->get_H0();
  const double OmegaR   = cosmo->get_OmegaR();  //OMEGA0 OR OMEGA(x)???
  const double OmegaCDM = cosmo->get_OmegaCDM();
  const double OmegaB   = cosmo->get_OmegaB();
  const double dtaudx   = rec->dtaudx_of_x(x);
  const double ddtauddx = rec->ddtauddx_of_x(x);
  const double R  = 4*OmegaR/(3*OmegaB*a);


  // SET: Scalar quantities (Phi, delta, v, ...)
  // ...
  // ...
  // ...
  // const double Theta2 = -4.0*ck/(9.0*Hp*dtaudx)*Theta[1];  // Theta[2] doesn't exist in tc array, so we set it here.

  const double Psi    = -Phi - 12.0*H0*H0/(ck*ck*a*a)*OmegaR*Theta[2];

  dPhidx        = Psi - ck*ck/(3.0*Hp*Hp)*Phi + H0*H0/(2.0*Hp*Hp)*(OmegaCDM/a*delta_cdm + OmegaB/a*delta_b + 4.0*OmegaR/(a*a)*Theta[0]);
  ddelta_cdmdx  = ck*v_cdm/(Hp) - 3.0*dPhidx;
  dv_cdmdx      = -v_cdm - ck/Hp*Psi;
  ddelta_bdx    = ck/Hp*v_b - 3.0*dPhidx;
  dv_bdx        = -v_b - ck/Hp*Psi + dtaudx*R*(3*Theta[1] + v_b);
  dThetadx[0]   = -ck/Hp*Theta[1] - dPhidx;
  dThetadx[1]   = ck/(3.0*Hp)*Theta[0] - 2.0*ck/(3.0*Hp)*Theta[2] + ck/(3.0*Hp)*Psi + dtaudx*(Theta[1] + 1.0/3.0*v_b);
  for(int l=2; l<Constants.n_ell_theta-1; l++){
    double kd = (l == 2) ? 1 : 0;  // kronecker delta for l==2.
    dThetadx[l] = l*ck/((2*l + 1)*Hp)*Theta[l-1] - (l + 1)*ck/((2*l + 1)*Hp)*Theta[l+1] + dtaudx*(Theta[l] - 0.1*Theta[2]*kd);
  }
  int l = Constants.n_ell_theta-1;
  dThetadx[l] = ck/Hp*Theta[l-1] - Constants.c*(l+1.0)/(Hp*cosmo->eta_of_x(x))*Theta[l] + dtaudx*Theta[l];


  // const double q = (-((1 - R)*dtaudx + (1 + R)*ddtauddx)*(3*Theta[1] + v_b) - ck/Hp*Psi + (1 - dHpdx/Hp)*ck/Hp*(-Theta[0] + 2*Theta2) - ck/Hp*dThetadx[0])/((1 + R)*dtaudx + dHpdx/Hp - 1);
  // dv_bdx = 1/(1 + R)*(-v_b - ck/Hp*Psi + R*(q + ck/Hp*(-Theta[0] + 2*Theta2) - ck/Hp*Psi));
  // dThetadx[1] = 1.0/3.0*(q - dv_bdx);


  return GSL_SUCCESS;
}

//====================================================
// Get methods
//====================================================

double Perturbations::get_delta_cdm(const double x, const double k) const{
  double log_k = log(k);
  return delta_cdm_spline(x,log_k);
}
double Perturbations::get_delta_b(const double x, const double k) const{
  double log_k = log(k);
  return delta_b_spline(x,log_k);
}
double Perturbations::get_v_cdm(const double x, const double k) const{
  double log_k = log(k);
  return v_cdm_spline(x,log_k);
}
double Perturbations::get_v_b(const double x, const double k) const{
  double log_k = log(k);
  return v_b_spline(x,log_k);
}
double Perturbations::get_Phi(const double x, const double k) const{
  double log_k = log(k);
  return Phi_spline(x,log_k);
}
double Perturbations::get_Psi(const double x, const double k) const{
  double log_k = log(k);
  return Psi_spline(x,log_k);
}
double Perturbations::get_Pi(const double x, const double k) const{
  double log_k = log(k);
  return Pi_spline(x,log_k);
}
double Perturbations::get_Source_T(const double x, const double k) const{
  double log_k = log(k);
  return ST_spline(x,log_k);
}
double Perturbations::get_Source_E(const double x, const double k) const{
  double log_k = log(k);
  return SE_spline(x,log_k);
}
double Perturbations::get_Theta(const double x, const double k, const int ell) const{
  double log_k = log(k);

  switch(ell){
    case 0:
      return Theta0_spline(x,log_k);
    case 1:
      return Theta1_spline(x,log_k);
    case 2:
      return Theta2_spline(x,log_k);
    default:
      printf("Error: Unknown element of theta, l=%d", ell);
      return -9999999;
  }

}


//====================================================
// Print some useful info about the class
//====================================================

void Perturbations::info() const{
  std::cout << "\n";
  std::cout << "Info about perturbations class:\n";
  std::cout << "x_start:       " << x_start                << "\n";
  std::cout << "x_end:         " << x_end                  << "\n";
  std::cout << "n_x:     " << n_x              << "\n";
  std::cout << "k_min (1/Mpc): " << k_min * Constants.Mpc  << "\n";
  std::cout << "k_max (1/Mpc): " << k_max * Constants.Mpc  << "\n";
  std::cout << "n_k:     " << n_k              << "\n";
  if(Constants.polarization)
    std::cout << "We include polarization\n";
  else
    std::cout << "We do not include polarization\n";
  if(Constants.neutrinos)
    std::cout << "We include neutrinos\n";
  else
    std::cout << "We do not include neutrinos\n";

  std::cout << "Information about the perturbation system:\n";
  std::cout << "ind_deltacdm:       " << Constants.ind_deltacdm         << "\n";
  std::cout << "ind_deltab:         " << Constants.ind_deltab           << "\n";
  std::cout << "ind_v_cdm:          " << Constants.ind_vcdm             << "\n";
  std::cout << "ind_v_b:            " << Constants.ind_vb               << "\n";
  std::cout << "ind_Phi:            " << Constants.ind_Phi              << "\n";
  std::cout << "ind_start_theta:    " << Constants.ind_start_theta      << "\n";
  std::cout << "n_ell_theta:        " << Constants.n_ell_theta          << "\n";
  if(Constants.polarization){
    std::cout << "ind_start_thetap:   " << Constants.ind_start_thetap   << "\n";
    std::cout << "n_ell_thetap:       " << Constants.n_ell_thetap       << "\n";
  }
  if(Constants.neutrinos){
    std::cout << "ind_start_nu:       " << Constants.ind_start_nu       << "\n";
    std::cout << "n_ell_neutrinos     " << Constants.n_ell_neutrinos    << "\n";
  }
  std::cout << "n_ell_tot_full:     " << Constants.n_ell_tot_full       << "\n";

  std::cout << "Information about the perturbation system in tight coupling:\n";
  std::cout << "ind_deltacdm:       " << Constants.ind_deltacdm_tc      << "\n";
  std::cout << "ind_deltab:         " << Constants.ind_deltab_tc        << "\n";
  std::cout << "ind_v_cdm:          " << Constants.ind_vcdm_tc          << "\n";
  std::cout << "ind_v_b:            " << Constants.ind_vb_tc            << "\n";
  std::cout << "ind_Phi:            " << Constants.ind_Phi_tc           << "\n";
  std::cout << "ind_start_theta:    " << Constants.ind_start_theta_tc   << "\n";
  std::cout << "n_ell_theta:        " << Constants.n_ell_theta_tc       << "\n";
  if(Constants.neutrinos){
    std::cout << "ind_start_nu:       " << Constants.ind_start_nu_tc    << "\n";
    std::cout << "n_ell_neutrinos     " << Constants.n_ell_neutrinos_tc << "\n";
  }
  std::cout << "n_ell_tot_tc:       " << Constants.n_ell_tot_tc         << "\n";
  std::cout << std::endl;
}

//====================================================
// Output some results to file for a given value of k
//====================================================

void Perturbations::output(const double k, const std::string filename) const{
  std::ofstream fp(filename.c_str());
  const int npts = 5000;
  auto x_array = Utils::linspace(x_start, 0, npts);
  auto print_data = [&] (const double x) {
    double arg = k * Constants.c * (cosmo->eta_of_x(0.0) - cosmo->eta_of_x(x));

    fp << std::setprecision(10);
    fp << x                  << " ";
    fp << get_delta_cdm(x,k) << " ";
    fp << get_delta_b(x,k)   << " ";
    fp << get_v_cdm(x,k)     << " ";
    fp << get_v_b(x,k)       << " ";
    fp << get_Theta(x,k,0)   << " ";
    fp << get_Theta(x,k,1)   << " ";
    fp << get_Theta(x,k,2)   << " ";
    fp << get_Phi(x,k)       << " ";
    fp << get_Psi(x,k)       << " ";
    fp << get_Pi(x,k)        << " ";
    // fp << get_Source_T(x,k)  << " ";
    // fp << get_Source_T(x,k) * Utils::j_ell(5,   arg)           << " ";
    // fp << get_Source_T(x,k) * Utils::j_ell(50,  arg)           << " ";
    // fp << get_Source_T(x,k) * Utils::j_ell(500, arg)           << " ";
    fp << "\n";
  };
  std::for_each(x_array.begin(), x_array.end(), print_data);
}

