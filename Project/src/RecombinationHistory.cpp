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
  Vector x_array = Utils::linspace(x_start, x_end, npts_rec_arrays);
  Vector Xe_arr(npts_rec_arrays);
  Vector ne_arr(npts_rec_arrays);

  // Calculate recombination history
  bool saha_regime = true;
  for(int i = 0; i < npts_rec_arrays; i++){

    //==============================================================
    // Get X_e from solving the Saha equation so
    // implement the function electron_fraction_from_saha_equation
    //==============================================================
    auto Xe_ne_data = electron_fraction_from_saha_equation(x_array[i]);

    // Electron fraction and number density
    const double Xe_current = Xe_ne_data.first;
    const double ne_current = Xe_ne_data.second;

    // Are we still in the Saha regime?
    if(Xe_current < Xe_saha_limit)
      saha_regime = false;

    if(saha_regime){
      
      //=============================================================================
      // Store the result we got from the Saha equation
      //=============================================================================
      Xe_arr[i] = Xe_current;
      ne_arr[i] = ne_current;
    } else {

      //==============================================================
      // Compute X_e from current time til today by solving the Peebles equation
      //==============================================================

      // The Peebles ODE equation
      ODESolver peebles_Xe_ode;
      ODEFunction dXedx = [&](double x, const double *Xe, double *dXedx){
        return rhs_peebles_ode(x, Xe, dXedx);
      };
      
      //=============================================================================
      // TODO: Set up IC, solve the ODE and fetch the result 
      //=============================================================================
      int remaining_indices = npts_rec_arrays - i;  // How many indices remains to be solved by Peebles?
      Vector x_array_Peebles(remaining_indices);
      for(int j=0; j<remaining_indices; j++){
        x_array_Peebles[j] = x_array[j+i];
      }
      Vector Xe_inc = {Xe_current};  // Initial condition of Peebles is last element of Saha solution.

      ODESolver ode(1e-8, 1e-12, 1e-12);
      ode.solve(dXedx, x_array_Peebles, Xe_inc);

      Vector Xe_Peebles = ode.get_data_by_component(0);

      const double OmegaB = cosmo->get_OmegaB();
      const double rho_c  = cosmo->get_rho_crit();
      const double m_H    = Constants.m_H;
      
      for(int j=0; j<remaining_indices; j++){
        Xe_arr[j+i] = Xe_Peebles[j];  // Putting the solutions from Peebles back into the Xe array.
        double a = exp(x_array_Peebles[j]);
        double n_H = OmegaB*rho_c/(m_H*a*a*a);
        ne_arr[j+i] = Xe_arr[j+i]*n_H;
      }

      break;
    }
  }

  //=============================================================================
  // Spline the result.
  //=============================================================================

  // Log(Xe) behaves more smoothly in x, so we spline the log, and simply exp it when we need it.
  Vector Xe_log_arr(npts_rec_arrays);
  for(int i=0; i<npts_rec_arrays; i++){
    if(Xe_arr[i] < exp(-20)){  // Setting a lower bound on Xe, as log(0)=-inf.
      Xe_log_arr[i] = -20;
    }
    else{
      Xe_log_arr[i] = log(Xe_arr[i]);
    }
  }

  log_Xe_of_x_spline.create(x_array, Xe_log_arr);

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
  const double OmegaB = cosmo->get_OmegaB();
  const double rho_c  = cosmo->get_rho_crit();
  const double T_b    = cosmo->get_TCMB()/a;  // Assuming T_b = T_CMB/a to be a good approximation.
  
  const double n_b    = OmegaB*rho_c/(m_H*a*a*a);
  const double n_H    = n_b;

  // Electron fraction and number density
  double Xe = 0.0;
  double ne = 0.0;
  
  //=============================================================================
  // Compute Xe and ne from the Saha equation
  //=============================================================================

  double A = 1/(n_b) * pow((m_e*k_b*T_b)/(2*M_PI*hbar*hbar), 1.5) * exp(-epsilon_0/(k_b*T_b));

  if(A < 1e-20){  // The solution to Xe becomes numerically unstable for large and small A,
    Xe = 0;       // so we hardcode the asymptotic solutions in each direction.
  }
  else if(A > 1e6){
    Xe = 1;
  }
  else{
    Xe = -A/2 + A/2*sqrt(1 + 4/A);
  }
  ne = Xe*n_H;

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
  const double OmegaB    = cosmo->get_OmegaB();
  const double T_b       = cosmo->get_TCMB()/a;
  const double rho_c     = cosmo->get_rho_crit();
  const double H         = cosmo->H_of_x(x);
  const double n_H = OmegaB*rho_c/(Constants.m_H*a*a*a);

  // Commonly used combination of units, predifed for speed and prevention of loss of numerical precision.
  const double c_hbar     = c*hbar;
  const double ep0_c_hbar = epsilon_0/c_hbar;
  const double kb_Tb      = k_b*T_b;
  const double ep0_kb_Tb  = epsilon_0/kb_Tb;

  // Expression in the Peebles RHS equation.
  const double phi_2         = 0.448*log(ep0_kb_Tb);
  const double n_1s          = (1 - X_e)*n_H;
  const double Lambda_alpha = H*27*ep0_c_hbar*ep0_c_hbar*ep0_c_hbar/(64*M_PI*M_PI*n_1s);
  const double Lambda_2s1s  = 8.227;
  const double alpha_2      = 8/sqrt(3*M_PI)*sigma_T*c*sqrt(ep0_kb_Tb)*phi_2;
  const double beta         = alpha_2*pow((m_e*kb_Tb/(2*M_PI*hbar*hbar)), 1.5)*exp(-ep0_kb_Tb);
  // expression for beta has been inserted right into beta_2, to diminish the effect of the large exponents.
  const double beta_2       = alpha_2*pow((m_e*kb_Tb/(2*M_PI*hbar*hbar)), 1.5)*exp(-0.25*ep0_kb_Tb);
  const double Cr           = (Lambda_2s1s + Lambda_alpha)/(Lambda_2s1s + Lambda_alpha + beta_2);

  dXedx[0] = Cr/H*(beta*(1 - X_e) - n_H*alpha_2*X_e*X_e);

  return GSL_SUCCESS;
}

//====================================================
// Solve for the optical depth tau, compute th visibility function and spline the result
//====================================================

void RecombinationHistory::solve_for_optical_depth_tau(){
  Utils::StartTiming("opticaldepth");

  const double c = Constants.c;
  const double sigma_T = Constants.sigma_T;

  // Set up x-arrays to integrate over.
  const int npts = npts_rec_arrays;
  const int x0_index = (int) (-x_start/((double) (x_end - x_start)/npts));  // The index closest to x=0.
  Vector x_array = Utils::linspace(x_start, 0, npts);
  // We in reality solve the ODE for tau in [0, -x_start], for increased numerical accuracy:
  Vector x_array_backwards = Utils::linspace(0, -x_start, npts);

  // The ODE system dtau/dx, dtau_noreion/dx and dtau_baryon/dx
  ODEFunction dtaudx = [&](double x, const double *tau, double *dtaudx){
    double H = cosmo->H_of_x(-x);
    double n_e = ne_of_x(-x);

    dtaudx[0] = -c*n_e*sigma_T/H;

    return GSL_SUCCESS;
  };

  //=============================================================================
  // Set up and solve the ODE and make tau splines
  //=============================================================================
  Vector tau_inc = {0};

  ODESolver ode(1e-8, 1e-12, 1e-12);
  ode.solve(dtaudx, x_array_backwards, tau_inc);
  Vector tau_arr_backwards = ode.get_data_by_component(0);
  Vector dtaudx_arr_backwards = ode.get_derivative_data_by_component(0);

  Vector tau_arr(npts);
  Vector dtaudx_arr(npts);
  for(int i=0; i<npts; i++){
    tau_arr[i] = -tau_arr_backwards[npts-i-1];  // Reverse direction of array, as we've solved it backwards.
    dtaudx_arr[i] = dtaudx_arr_backwards[npts-i-1];
  }
  
  tau_of_x_spline.create(x_array, tau_arr);
  dtaudx_of_x_spline.create(x_array, dtaudx_arr);
  //=============================================================================
  // Compute visibility functions and spline everything
  //=============================================================================
  Vector g_arr(npts);
  double g_integral = 0;
  for(int i=0; i<npts; i++){
    g_arr[i] = -dtaudx_arr[i]*exp(-tau_arr[i]);
    g_integral += g_arr[i];
  }
  g_tilde_of_x_spline.create(x_array, g_arr);
  g_integral *= (0 - x_start)/npts;
  printf("Integral over g_tilde = %f\n", g_integral);

  Utils::EndTiming("opticaldepth");
}

//====================================================
// Get methods
//====================================================

double RecombinationHistory::tau_of_x(double x) const{
  // The spline only exists up until x=0, but tau is, by definition, 0 after that.
  if(x < 0){
    return tau_of_x_spline(x);
  }
  else{
    return 0.0;
  }
}

double RecombinationHistory::dtaudx_of_x(double x) const{
  // The spline only exists up until x=0, but dtaudx is, by definition, 0 after that.
  if(x < 0){
    return dtaudx_of_x_spline(x);
  }
  else{
    return 0.0;
  }
}

double RecombinationHistory::ddtauddx_of_x(double x) const{
  // The spline only exists up until x=0, but ddtauddx is, by definition, 0 after that.
  if(x < 0){
    return dtaudx_of_x_spline.deriv_x(x);
  }
  else{
    return 0.0;
  }
}

double RecombinationHistory::g_tilde_of_x(double x) const{
  // The spline only exists up until x=0, but g_tilde is, by definition, 0 after that.
  if(x < 0){
    return g_tilde_of_x_spline(x);
  }
  else{
    return 0.0;
  }
}

double RecombinationHistory::dgdx_tilde_of_x(double x) const{
  // The spline only exists up until x=0, but dgdx_tilde is, by definition, 0 after that.
  if(x < 0){
    return g_tilde_of_x_spline.deriv_x(x);
  }
  else{
    return 0.0;
  }
}

double RecombinationHistory::ddgddx_tilde_of_x(double x) const{
  // The spline only exists up until x=0, but ddgddx_tilde is, by definition, 0 after that.
  if(x < 0){
    return g_tilde_of_x_spline.deriv_xx(x);
  }
  else{
    return 0.0;
  }
}

double RecombinationHistory::Xe_of_x(double x) const{
  // The spline is stored as a log, so we exp it to get Xe.
  double log_Xe_of_x = log_Xe_of_x_spline(x);
  return exp(log_Xe_of_x);
}

double RecombinationHistory::ne_of_x(double x) const{
  double OmegaB = cosmo->get_OmegaB();
  double rho_c  = cosmo->get_rho_crit();
  double a = exp(x);
  const double n_H = OmegaB*rho_c/(Constants.m_H*a*a*a);
  return Xe_of_x(x)*n_H;
}

double RecombinationHistory::get_Yp() const{
  return Yp;
}

void RecombinationHistory::set_saha_limit(double saha_limit){
  // Function for setting the Saha limit to something else than the default 0.99.
  Xe_saha_limit = saha_limit;
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
  const int npts       = 100000;
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

