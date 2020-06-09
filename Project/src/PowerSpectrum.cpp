#include"PowerSpectrum.h"

//====================================================
// Constructors
//====================================================

PowerSpectrum::PowerSpectrum(
    BackgroundCosmology *cosmo, 
    RecombinationHistory *rec, 
    Perturbations *pert) : 
  cosmo(cosmo), 
  rec(rec), 
  pert(pert)
{}

//====================================================
// Do all the solving
//====================================================
void PowerSpectrum::solve(){

  //=========================================================================
  // TODO: Choose the range of k's and the resolution to compute Theta_ell(k)
  //=========================================================================
  Vector log_k_array = Utils::linspace(log(k_min), log(k_max), n_k);
  Vector k_array = exp(log_k_array);

  //=========================================================================
  // TODO: Make splines for j_ell. 
  // Implement generate_bessel_function_splines
  //=========================================================================
  generate_bessel_function_splines();

  //=========================================================================
  // TODO: Line of sight integration to get Theta_ell(k)
  // Implement line_of_sight_integration
  //=========================================================================
  line_of_sight_integration(k_array);

  //=========================================================================
  // TODO: Integration to get Cell by solving dCell^f/dlogk = Delta(k) * f_ell(k)^2
  // Implement solve_for_cell
  //=========================================================================
  auto cell_TT = solve_for_cell(log_k_array, thetaT_ell_of_k_spline, thetaT_ell_of_k_spline);
  cell_TT_spline.create(ells, cell_TT, "Cell_TT_of_ell");
  
  //=========================================================================
  // TODO: Do the same for polarization...
  //=========================================================================
  // ...
  // ...
  // ...
  // ...
}

//====================================================
// Generate splines of j_ell(z) needed for LOS integration
//====================================================

void PowerSpectrum::generate_bessel_function_splines(){
  Utils::StartTiming("besselspline");
  
  // Make storage for the splines
  j_ell_splines = std::vector<Spline>(ells.size());
    
  //=============================================================================
  // TODO: Compute splines for bessel functions j_ell(z)
  // Choose a suitable range for each ell
  // NB: you don't want to go larger than z ~ 40000, then the bessel routines
  // might break down. Use j_ell(z) = Utils::j_ell(ell, z)
  //=============================================================================

  for(size_t i = 0; i < ells.size(); i++){
    const int ell = ells[i];

    Vector log_z_array = Utils::linspace(log(1e-8), log(4e4), 10001);
    Vector z_array = exp(log_z_array);
    // Vector z_array = Utils::linspace(0, 4e4, 100001);
    Vector j_ell_array(10001);
    for(int j=0; j<10001; j++){
      j_ell_array[j] = Utils::j_ell(ell, z_array[j]);
    }

    j_ell_splines[i] = Spline(log_z_array, j_ell_array, "j_ell_spline");

  }

  Utils::EndTiming("besselspline");
}

//====================================================
// Do the line of sight integration for a single
// source function
//====================================================

Vector2D PowerSpectrum::line_of_sight_integration_single(
    Vector & k_array, 
    std::function<double(double,double)> &source_function){
  Utils::StartTiming("lineofsight");

  // Make storage for the results
  Vector2D result = Vector2D(ells.size(), Vector(k_array.size()));
  int nx = 2001;

  #pragma omp parallel for schedule(dynamic, 1)
  for(size_t ik = 0; ik < k_array.size(); ik++){

    //=============================================================================
    // TODO: Implement to solve for the general line of sight integral 
    // F_ell(k) = Int dx jell(k(eta-eta0)) * S(x,k) for all the ell values for the 
    // given value of k
    //=============================================================================
    // ...
    // ...
    // ...
    double k = k_array[ik];
    Vector x_array = Utils::linspace(-16, 0, nx);
    Vector Theta_ell_ini = {0.0};


    // ---- ODE ----
    Spline j_ell_spline;
    for(int iell = 0; iell < ells.size(); iell++){
      int ell = ells[iell];
      j_ell_spline = j_ell_splines[iell];


      // ---- integral ----
      double integral = 0;
      for(int ix=0; ix<nx; ix++){
        double x = x_array[ix];
        double Delta_eta = cosmo->eta_of_x(0) - cosmo->eta_of_x(x);
        double S = pert->get_Source_T(x,k);
        double term = k*Delta_eta;
        if(term < 1e-8){
          term = 1e-8;
        }
        double j = j_ell_spline(log(term));
        integral += j*S;
      }
      integral *= 16.0/nx;
      result[iell][ik] = integral;

    }
    // Store the result for Source_ell(k) in results[ell][ik]
  }

  // myfile.close();
  Utils::EndTiming("lineofsight");
  return result;
}

//====================================================
// Do the line of sight integration
//====================================================
void PowerSpectrum::line_of_sight_integration(Vector & k_array){
  const int n_k        = k_array.size();
  const int nells      = ells.size();
  
  // Make storage for the splines we are to create
  thetaT_ell_of_k_spline = std::vector<Spline>(nells);

  //============================================================================
  // TODO: Solve for Theta_ell(k) and spline the result
  //============================================================================

  // Make a function returning the source function
  std::function<double(double,double)> source_function_T = [&](double x, double k){
    return pert->get_Source_T(x,k);
  };

  // Do the line of sight integration
  Vector2D thetaT_ell_of_k = line_of_sight_integration_single(k_array, source_function_T);

  // Spline the result and store it in thetaT_ell_of_k_spline
  // ...
  // ...
  // ...
  // ...
  for(int iell=0; iell<ells.size(); iell++){
    Spline thetaT_spline = Spline(k_array, thetaT_ell_of_k[iell], "theta_T_spline");
    thetaT_ell_of_k_spline[iell] = thetaT_spline;
  }



}

//====================================================
// Compute Cell (could be TT or TE or EE) 
// Cell = Int_0^inf 4 * pi * P(k) f_ell g_ell dk/k
//====================================================
Vector PowerSpectrum::solve_for_cell(
    Vector & log_k_array,
    std::vector<Spline> & f_ell_spline,
    std::vector<Spline> & g_ell_spline){
  const int nells      = ells.size();

  //============================================================================
  // TODO: Integrate Cell = Int 4 * pi * P(k) f_ell g_ell dk/k
  // or equivalently solve the ODE system dCell/dlogk = 4 * pi * P(k) * f_ell * g_ell
  //============================================================================
  Vector result(nells);

  for(int iell=0; iell<nells; iell++){
    ODEFunction dCelldlogk_func = [&](double logk, const double *Cell, double *dCelldlogk){
      const double k = exp(logk);
      const double f_ell = f_ell_spline[iell](k);
      const double g_ell = g_ell_spline[iell](k);
      dCelldlogk[0] = 4*M_PI*A_s*pow(k/kpivot, n_s-1)*f_ell*g_ell;

      return GSL_SUCCESS;
    };

    Vector logk_array = Utils::linspace(log(k_min), log(k_max), n_k);
    Vector Cell_ic = {0.0};

    ODESolver ode;
    ode.solve(dCelldlogk_func, logk_array, Cell_ic);

    result[iell] = ode.get_data_by_component(0)[n_k-1];
  }

  return result;
}


//====================================================
// Get methods
//====================================================
double PowerSpectrum::get_cell_TT(const double ell) const{
  return cell_TT_spline(ell);
}
double PowerSpectrum::get_cell_TE(const double ell) const{
  return cell_TE_spline(ell);
}
double PowerSpectrum::get_cell_EE(const double ell) const{
  return cell_EE_spline(ell);
}
double PowerSpectrum::get_P(const double x, const double k) const{
  const double a = exp(x);
  const double ck = Constants.c*k;
  const double H0 = cosmo->get_H0();
  const double OmegaM = cosmo->get_OmegaB(x) + cosmo->get_OmegaCDM(x);
  const double Delta_m = ck*ck*pert->get_Psi(x, k)/(1.5*OmegaM/a*H0*H0);
  const double P_primordial = 2*M_PI*M_PI/(k*k*k)*A_s*pow(k/kpivot, n_s-1);
  return Delta_m*Delta_m*P_primordial;
}
double PowerSpectrum::get_Theta_ell_k(const double iell, const double k) const{
  return thetaT_ell_of_k_spline[iell](k);
}

//====================================================
// Output the cells to file
//====================================================

void PowerSpectrum::output_P(std::string filename) const{
  std::ofstream fp(filename.c_str());
  Vector logk_array = Utils::linspace(log(k_min), log(k_max), n_k);
  auto print_data = [&] (const double logk) {
    const double k = exp(logk);
    fp << k << " ";
    fp << get_P(0, k);
    fp << "\n";
  };
  std::for_each(logk_array.begin(), logk_array.end(), print_data);
}

void PowerSpectrum::output_Theta(std::string filename, const int iell) const{
  std::ofstream fp(filename.c_str());
  Vector logk_array = Utils::linspace(log(k_min), log(k_max), n_k);
  auto print_data = [&] (const double logk) {
    const double k = exp(logk);
    fp << k << " ";
    fp << get_Theta_ell_k(iell, k);
    fp << "\n";
  };
  std::for_each(logk_array.begin(), logk_array.end(), print_data);
}

void PowerSpectrum::output_Cell(std::string filename) const{
  std::ofstream fp(filename.c_str());
  const int ellmax = int(ells[ells.size()-1]);
  auto ellvalues = Utils::linspace(2, ellmax, ellmax-1);
  auto print_data = [&] (const double ell) {
    double normfactor  = (ell * (ell+1)) / (2.0 * M_PI) * pow(1e6 * cosmo->get_TCMB(), 2);
    double normfactorN = (ell * (ell+1)) / (2.0 * M_PI) 
      * pow(1e6 * cosmo->get_TCMB() *  pow(4.0/11.0, 1.0/3.0), 2);
    double normfactorL = (ell * (ell+1)) * (ell * (ell+1)) / (2.0 * M_PI);
    fp << ell                                 << " ";
    fp << cell_TT_spline( ell ) * normfactor  << " ";
    if(Constants.polarization){
      fp << cell_EE_spline( ell ) * normfactor  << " ";
      fp << cell_TE_spline( ell ) * normfactor  << " ";
    }
    fp << "\n";
  };
  std::for_each(ellvalues.begin(), ellvalues.end(), print_data);
}

