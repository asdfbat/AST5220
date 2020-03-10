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
  const int x_start = -16;
  const int x_end = 6;
  Vector x_array = Utils::linspace(x_start, x_end, npts_rec_arrays);
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
      saha_regime = false;

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
      int remaining_indices = npts_rec_arrays - i;
      Vector x_array_Peebles(remaining_indices);
      for(int j=0; j<remaining_indices; j++){
        x_array_Peebles[j] = x_array[j+i];
      }
      Vector Xe_inc = {Xe_current};

      ODESolver ode;
      ode.solve(dXedx, x_array_Peebles, Xe_inc);

      Vector Xe_Peebles = ode.get_data_by_component(0);
      // for(int j=0; j<remaining_indices; j++){
      //   printf("%f %e\n", x_array[i+j], Xe_Peebles[j]);
      // }

      const double OmegaB = cosmo->get_OmegaB();
      const double rho_c  = cosmo->get_rho_crit();
      const double m_H    = Constants.m_H;
      
      for(int j=0; j<remaining_indices; j++){
        Xe_arr[j+i] = Xe_Peebles[j];
        double a = exp(x_array_Peebles[j]);
        double n_H = OmegaB*rho_c/(m_H*a*a*a);
        ne_arr[j+i] = Xe_arr[j]*n_H;
      }

      break;
    }
  }

  //=============================================================================
  // TODO: Spline the result. Implement and make sure the Xe_of_x, ne_of_x 
  // functions are working
  //=============================================================================
  //...
  //...
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
  // tau_of_x_spline.create(x_array, x_array); //PLACEHODLER: WRONG!
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

  if(A < 1e-30){  // The solution to Xe becomes numerically unstable for large and small A,
    Xe = 0;       // so we hardcode the asymptotic solutions in each direction.
  }
  else if(A > 1e12){
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
  const long double k_b         = Constants.k_b;
  const long double G           = Constants.G;
  const long double c           = Constants.c;
  const long double m_e         = Constants.m_e;
  const long double hbar        = Constants.hbar;
  const long double m_H         = Constants.m_H;
  const long double sigma_T     = Constants.sigma_T;
  const long double lambda_2s1s = Constants.lambda_2s1s;
  const long double epsilon_0   = Constants.epsilon_0;

  // Cosmological parameters
  // const double OmegaB      = cosmo->get_OmegaB();
  // ...
  // ...
  const long double OmegaB    = cosmo->get_OmegaB();
  const long double T_b       = cosmo->get_TCMB()/a;
  const long double rho_c     = cosmo->get_rho_crit();
  const long double H         = cosmo->H_of_x(x);
  const long double n_H = OmegaB*rho_c/(Constants.m_H*a*a*a);

  // Commonly used combination of units, predifed for speed and prevention of loss of numerical precision.
  const long double c_hbar     = c*hbar;
  const long double ep0_c_hbar = epsilon_0/c_hbar;
  const long double kb_Tb      = k_b*T_b;
  const long double ep0_kb_Tb  = epsilon_0/kb_Tb;

  // Expression in the Peebles RHS equation.
  const long double phi_2         = 0.448*log(ep0_kb_Tb);
  const long double n_1s          = (1 - X_e)*n_H;
  const long double Lambda_alpha = H*27*ep0_c_hbar*ep0_c_hbar*ep0_c_hbar/(64*M_PI*M_PI*n_1s);
  const long double Lambda_2s1s  = 8.227;
  const long double alpha_2      = 8/sqrt(3*M_PI)*sigma_T*c*sqrt(ep0_kb_Tb)*phi_2;
  const long double beta         = alpha_2*pow((m_e*kb_Tb/(2*M_PI*hbar*hbar)), 1.5)*exp(-ep0_kb_Tb);
  const long double beta_2       = alpha_2*pow((m_e*kb_Tb/(2*M_PI*hbar*hbar)), 1.5)*exp(-0.25*ep0_kb_Tb);
  const long double beta_22      = beta*exp(0.75*ep0_kb_Tb);
  const long double Cr           = (Lambda_2s1s + Lambda_alpha)/(Lambda_2s1s + Lambda_alpha + beta_2);

  //     hbar            c           ep_0          k_b         c*hbar       ep0_c_hbar
  // 1.054572e-34  2.997925e+08  2.179872e-18  1.380649e-23  3.161527e-26  6.894998e+07
  
  //           kB_Tb                     ep0_kB_Tb
  // [6.0900e-20 - 9.3257e-26]   [3.5794e+01 - 2.3375e+07]


  // static bool asdf = true;
  // if(asdf){
  //   printf("%10s  %10s  %10s  %10s  %10s  %10s  %10s  %10s\n", "x", "Tb", "alpha2", "ep0_kb_Tb", "phi2", "beta", "beta_2", "beta_22");
  //   asdf=false;
  // }

  // printf("%10.4e  %10.4e  %10.4e  %10.4e  %10.4e  %10.4e,  %10.4e  %10.4e  %10.4e  %10.4e\n", x, (double) ep0_kb_Tb, (double) alpha_2, (double) ep0_kb_Tb, (double) phi_2, (double) beta, (double) beta_2, (double) beta_22, (double) alpha_2*pow((m_e*kb_Tb/(2*M_PI*hbar*hbar)), 1.5), (double) exp(-1/4*ep0_kb_Tb));

  //=============================================================================
  // TODO: Write the expression for dXedx
  //=============================================================================
  //...
  //...
  // if(X_e > 0.98){
  //   printf("%f %e %e\n", x, X_e, Cr/H*(beta*(1 - X_e) - n_H*alpha_2*X_e*X_e));
  // }
  dXedx[0] = Cr/H*(beta*(1 - X_e) - n_H*alpha_2*X_e*X_e);

  // printf("%e %e %e %e %e %e %e %e\n", X_e, Cr, H, n_H, alpha_2, beta, beta_2, phi_2);
  // printf("%e %e %e %e %e %e\n", x, Xe[0], dXedx[0], Cr/H, beta*(1 - X_e), n_H*alpha_2*X_e*X_e);
  return GSL_SUCCESS;
}

//====================================================
// Solve for the optical depth tau, compute the 
// visibility function and spline the result
//====================================================

void RecombinationHistory::solve_for_optical_depth_tau(){
  Utils::StartTiming("opticaldepth");

  const double c = Constants.c;
  const double sigma_T = Constants.sigma_T;

  // Set up x-arrays to integrate over. We split into three regions as we need extra points in reionisation
  const int npts = 1e5;
  const int x_start = -16;
  const int x_end = 4;
  const int x0_index = (int) (-x_start/((double) (x_end - x_start)/npts));
  Vector x_array = Utils::linspace(x_start, x_end, npts);

  // The ODE system dtau/dx, dtau_noreion/dx and dtau_baryon/dx
  ODEFunction dtaudx = [&](double x, const double *tau, double *dtaudx){

    //=============================================================================
    // TODO: Write the expression for dtaudx
    //=============================================================================
    //...
    //...
    double H = cosmo->H_of_x(x);
    double n_e = ne_of_x(x);

    dtaudx[0] = -c*n_e*sigma_T/H;

    return GSL_SUCCESS;
  };

  //=============================================================================
  // TODO: Set up and solve the ODE and make tau splines
  //=============================================================================
  //...
  //...

  Vector tau_inc = {1e5};

  ODESolver ode;
  ode.solve(dtaudx, x_array, tau_inc);
  Vector tau_arr = ode.get_data_by_component(0);
  Vector dtaudx_arr = ode.get_derivative_data_by_component(0);


  double tau_0 = tau_arr[x0_index];
  for(int i=0; i<npts; i++){
    tau_arr[i] -= tau_0;
    if(i%100 == 0)
    printf("%e\n", dtaudx_arr[i]);
  }

  tau_of_x_spline.create(x_array, tau_arr);
  dtaudx_of_x_spline.create(x_array, dtaudx_arr);
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

  return dtaudx_of_x_spline(x);
}

double RecombinationHistory::ddtauddx_of_x(double x) const{

  //=============================================================================
  // TODO: Implement
  //=============================================================================
  //...
  //...

  return dtaudx_of_x_spline.deriv_x(x);
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
  const int npts       = 5000;
  const double x_min   = -16;
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

