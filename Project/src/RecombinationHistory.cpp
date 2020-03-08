#include"RecombinationHistory.h"

//====================================================
// Constructors
//====================================================
   
RecombinationHistory::RecombinationHistory(
    BackgroundCosmology *cosmo, 
    double Yp) :
  cosmo(cosmo),
  Yp(Yp)
{}

//====================================================
// Do all the solving we need to do
//====================================================

void RecombinationHistory::solve(){
    
  // Compute and spline Xe, ne
  solve_number_density_electrons();
   
  // Compute and spline tau, dtaudx, ddtauddx, g, dgdx, ddgddx, ...
  solve_for_optical_depth_tau();
}

//====================================================
// Solve for X_e and n_e using Saha and Peebles and spline the result
//====================================================

void RecombinationHistory::solve_number_density_electrons(){
  Utils::StartTiming("Xe");
  
  //=============================================================================
  // TODO: Set up x-array and make arrays to store X_e(x) and n_e(x) on
  //=============================================================================
  Vector x_array = Utils::linspace(Constants.x_start, Constants.x_end, npts_rec_arrays);
  Vector Xe_arr(npts_rec_arrays);
  Vector ne_arr(npts_rec_arrays);

  // Calculate recombination history
  bool saha_regime = true;
  for(int i = 0; i < npts_rec_arrays; i++){

    //==============================================================
    // TODO: Get X_e from solving the Saha equation so
    // implement the function electron_fraction_from_saha_equation
    //==============================================================
    auto Xe_ne_data = electron_fraction_from_saha_equation(x_array[i]);

    // Electron fraction and number density
    const double Xe_current = Xe_ne_data.first;
    const double ne_current = Xe_ne_data.second;

    // Are we still in the Saha regime?
    if(Xe_current < Xe_saha_limit)
      // saha_regime = false;
      saha_regime = true;

    if(saha_regime){
      
      //=============================================================================
      // TODO: Store the result we got from the Saha equation
      //=============================================================================
      //...
      //...
      Xe_arr[i] = Xe_current;
      ne_arr[i] = ne_current;

    } else {

      //==============================================================
      // TODO: Compute X_e from current time til today by solving 
      // the Peebles equation (NB: if you solve all in one go remember to
      // exit the for-loop!)
      // Implement rhs_peebles_ode
      //==============================================================
      //...
      //...

      // The Peebles ODE equation
      ODESolver peebles_Xe_ode;
      ODEFunction dXedx = [&](double x, const double *Xe, double *dXedx){
        return rhs_peebles_ode(x, Xe, dXedx);
      };
      
      //=============================================================================
      // TODO: Set up IC, solve the ODE and fetch the result 
      //=============================================================================
      //...
      //...
    
    }
  }

  //=============================================================================
  // TODO: Spline the result. Implement and make sure the Xe_of_x, ne_of_x 
  // functions are working
  //=============================================================================
  //...
  //...
  Vector Xe_log_arr(npts_rec_arrays);
  for(int i=0; i<npts_rec_arrays; i++){
    if(Xe_arr[i] < exp(-16)){
      Xe_log_arr[i] = -16;
    }
    else{
      Xe_log_arr[i] = log(Xe_arr[i]);
    }
    // printf("%f\n", Xe_log_arr[i]);
  }

  log_Xe_of_x_spline.create(x_array, Xe_log_arr);
  tau_of_x_spline.create(x_array, x_array); //PLACEHODLER: WRONG!
  g_tilde_of_x_spline.create(x_array, x_array); //PLACEHOLDER: WRONG!

  Utils::EndTiming("Xe");
}

//====================================================
// Solve the Saha equation to get ne and Xe
//====================================================
std::pair<double,double> RecombinationHistory::electron_fraction_from_saha_equation(double x) const{
  const double a           = exp(x);
 
  // Physical constants
  const double k_b         = Constants.k_b;
  const double G           = Constants.G;
  const double m_e         = Constants.m_e;
  const double hbar        = Constants.hbar;
  const double m_H         = Constants.m_H;
  const double epsilon_0   = Constants.epsilon_0;
  const double H0_over_h   = Constants.H0_over_h;

  // Fetch cosmological parameters
  //const double OmegaB      = cosmo->get_OmegaB();
  //...
  //...
  const double OmegaB = cosmo->get_OmegaB(x);
  const double rho_c  = cosmo->get_rho_crit();
  const double T_b    = cosmo->get_TCMB()/a;  // Assuming T_b = T_CMB/a to be a good approximation.
  
  const double n_b    = OmegaB*rho_c/(m_H*a*a*a);
  const double n_H    = n_b;

  // Electron fraction and number density
  double Xe = 0.0;
  double ne = 0.0;
  
  //=============================================================================
  // TODO: Compute Xe and ne from the Saha equation
  //=============================================================================
  //...
  //...

  // double log_n_b  = log(n_b);
  // double log_m_e  = log(m_e);
  // double log_k_b  = log(k_b);
  // double log_T_b  = log(T_b);
  // double log_hbar = log(hbar);

  // double log_A = 


  double A = 1/(n_b) * pow((m_e*k_b*T_b)/(2*M_PI*hbar*hbar), 1.5) * exp(-epsilon_0/(k_b*T_b));

  if(A < 1e-12){
    Xe = 0;
  }
  else if(A > 1e8){
    Xe = 1;
  }
  else{
    Xe = -A/2 + A/2*sqrt(1 + 4/A);
  }
  ne = Xe*n_H;

  // printf("%e %e\n", A, Xe);

  return std::pair<double,double>(Xe, ne);
}

//====================================================
// The right hand side of the dXedx Peebles ODE
//====================================================
int RecombinationHistory::rhs_peebles_ode(double x, const double *Xe, double *dXedx){

  // Current value of a and X_e
  const double X_e         = Xe[0];
  const double a           = exp(x);

  // Physical constants in SI units
  const double k_b         = Constants.k_b;
  const double G           = Constants.G;
  const double c           = Constants.c;
  const double m_e         = Constants.m_e;
  const double hbar        = Constants.hbar;
  const double m_H         = Constants.m_H;
  const double sigma_T     = Constants.sigma_T;
  const double lambda_2s1s = Constants.lambda_2s1s;
  const double epsilon_0   = Constants.epsilon_0;

  // Cosmological parameters
  // const double OmegaB      = cosmo->get_OmegaB();
  // ...
  // ...

  //=============================================================================
  // TODO: Write the expression for dXedx
  //=============================================================================
  //...
  //...
  
  dXedx[0] = 0.0;

  return GSL_SUCCESS;
}

//====================================================
// Solve for the optical depth tau, compute the 
// visibility function and spline the result
//====================================================

void RecombinationHistory::solve_for_optical_depth_tau(){
  Utils::StartTiming("opticaldepth");

  // Set up x-arrays to integrate over. We split into three regions as we need extra points in reionisation
  const int npts = 1000;
  Vector x_array = Utils::linspace(x_start, x_end, npts);

  // The ODE system dtau/dx, dtau_noreion/dx and dtau_baryon/dx
  ODEFunction dtaudx = [&](double x, const double *tau, double *dtaudx){

    //=============================================================================
    // TODO: Write the expression for dtaudx
    //=============================================================================
    //...
    //...

    // Set the derivative for photon optical depth
    dtaudx[0] = 0.0;

    return GSL_SUCCESS;
  };

  //=============================================================================
  // TODO: Set up and solve the ODE and make tau splines
  //=============================================================================
  //...
  //...

  //=============================================================================
  // TODO: Compute visibility functions and spline everything
  //=============================================================================
  //...
  //...

  Utils::EndTiming("opticaldepth");
}

//====================================================
// Get methods
//====================================================

double RecombinationHistory::tau_of_x(double x) const{
  return tau_of_x_spline(x);
}

double RecombinationHistory::dtaudx_of_x(double x) const{

  //=============================================================================
  // TODO: Implement. Either from the tau-spline tau_of_x_spline.deriv_(x) or 
  // from a separate spline if you choose to do this
  //=============================================================================
  //...
  //...

  return 0.0;
}

double RecombinationHistory::ddtauddx_of_x(double x) const{

  //=============================================================================
  // TODO: Implement
  //=============================================================================
  //...
  //...

  return 0.0;
}

double RecombinationHistory::g_tilde_of_x(double x) const{
  return g_tilde_of_x_spline(x);
}

double RecombinationHistory::dgdx_tilde_of_x(double x) const{

  //=============================================================================
  // TODO: Implement
  //=============================================================================
  //...
  //...

  return 0.0;
}

double RecombinationHistory::ddgddx_tilde_of_x(double x) const{

  //=============================================================================
  // TODO: Implement
  //=============================================================================
  //...
  //...

  return 0.0;
}

double RecombinationHistory::Xe_of_x(double x) const{

  //=============================================================================
  // TODO: Implement
  //=============================================================================
  //...
  //...
  double log_Xe_of_x = log_Xe_of_x_spline(x);
  return exp(log_Xe_of_x);
}

double RecombinationHistory::ne_of_x(double x) const{

  //=============================================================================
  // TODO: Implement
  //=============================================================================
  //...
  //...
  //...
  double OmegaB = cosmo->get_OmegaB(x);
  double rho_c  = cosmo->get_rho_crit();
  double a = exp(x);
  const double n_H = OmegaB*rho_c/(Constants.m_H*a*a*a);
  return Xe_of_x(x)*n_H;
}

double RecombinationHistory::get_Yp() const{
  return Yp;
}

//====================================================
// Print some useful info about the class
//====================================================
void RecombinationHistory::info() const{
  std::cout << "\n";
  std::cout << "Info about recombination/reionization history class:\n";
  std::cout << "Yp:          " << Yp          << "\n";
  std::cout << std::endl;
} 

//====================================================
// Output the data computed to file
//====================================================
void RecombinationHistory::output(const std::string filename) const{
  std::ofstream fp(filename.c_str());
  const int npts       = 5000;
  const double x_min   = x_start;
  const double x_max   = x_end;

  Vector x_array = Utils::linspace(x_min, x_max, npts);
  auto print_data = [&] (const double x) {
    fp << x                    << " ";  // 0
    fp << Xe_of_x(x)           << " ";  // 1
    fp << ne_of_x(x)           << " ";  // 2
    fp << tau_of_x(x)          << " ";  // 3
    fp << dtaudx_of_x(x)       << " ";  // 4
    fp << ddtauddx_of_x(x)     << " ";  // 5
    fp << g_tilde_of_x(x)      << " ";  // 6
    fp << dgdx_tilde_of_x(x)   << " ";  // 7
    fp << ddgddx_tilde_of_x(x) << " ";  // 8
    fp << "\n";
  };
  std::for_each(x_array.begin(), x_array.end(), print_data);
}

