#include "BackgroundCosmology.h"

//====================================================
// Constructors
//====================================================
    
BackgroundCosmology::BackgroundCosmology(
    double h, 
    double OmegaB, 
    double OmegaCDM, 
    double Neff,
    double TCMB) :
  h(h),
  OmegaB(OmegaB),
  OmegaCDM(OmegaCDM),
  Neff(Neff), 
  TCMB(TCMB)
{
  H0 = Constants.H0_over_h*h;
  rho_c = 3*H0*H0/(8*M_PI*Constants.G);
  OmegaR = pow(M_PI, 2)/15.0*pow(Constants.k_b*TCMB, 4)/(pow(Constants.hbar, 3)*pow(Constants.c, 5))*8*M_PI*Constants.G/(3*H0*H0);
  OmegaLambda = 1 - OmegaR - OmegaB - OmegaCDM;
  OmegaNu = 0;
  OmegaK = 0;
  printf("Sum of Omegas = %e\n", OmegaR + OmegaCDM + OmegaB + OmegaLambda);
  printf("OmegaR0 = %e\n", OmegaR);
  printf("H0 in SI: %e\n", H0);
}

// Solve the background
void BackgroundCosmology::solve(){
  Utils::StartTiming("Eta");
    
  int npts = 100000;
  double a_start = 1e-12;
  double a_stop = 1e2;
  double x_start = log(a_start);
  double x_stop = log(a_stop);
  printf("Solving conformal time for x in [%f, %f]\n", x_start, x_stop);

  Vector x_array = Utils::linspace(x_start, x_stop, npts);


  // The ODE for deta/dx
  ODEFunction detadx = [&](double x, const double *eta, double *detadx){
    detadx[0] = Constants.c/Hp_of_x(x);
    return GSL_SUCCESS;
  };

  //=============================================================================
  // Set the initial condition, set up the ODE system, solve and make the spline eta_of_x_spline.
  //=============================================================================

  Vector eta_inc = {0.0};

  ODESolver ode;
  ode.solve(detadx, x_array, eta_inc);
  Vector2D all_data = ode.get_data();
  Vector eta_data(all_data.size());
  for(int i=0; i<all_data.size(); i++){
    eta_data[i] = all_data[i][0];
  }

  eta_of_x_spline.create(x_array, eta_data, "eta");
  printf("%e %e \n", x_array[x_array.size()-1], eta_data[x_array.size()-1]);

  Utils::EndTiming("Eta");
}


//====================================================
// Get methods
//====================================================


double BackgroundCosmology::H_of_x(double x) const{  // Note: SI Units are used, not km/Mpc/s.
  double Omega_sum = OmegaCDM*exp(-3*x) + OmegaB*exp(-3*x) + OmegaR*exp(-4*x) + OmegaLambda;
  return get_H0()*sqrt(Omega_sum);
}


double BackgroundCosmology::Hp_of_x(double x) const{
  return exp(x)*H_of_x(x);
}


double BackgroundCosmology::dHpdx_of_x(double x) const{
  double exp_m1 = exp(-x);
  double exp_m2 = exp(-2*x);
  double exp_p2 = exp(2*x);

  double Omega_M = OmegaCDM + OmegaB;
  double Omega_sum1 = -Omega_M*exp_m1 - 2*OmegaR*exp_m2 + 2*OmegaLambda*exp_p2;
  double Omega_sum2 = Omega_M*exp_m1 + OmegaR*exp_m2 + OmegaLambda*exp_p2;
  return get_H0()*Omega_sum1/(2*sqrt(Omega_sum2));
}


double BackgroundCosmology::ddHpddx_of_x(double x) const{
  double exp_m1 = exp(-x);
  double exp_m2 = exp(-2*x);
  double exp_p2 = exp(2*x);

  double Omega_M = OmegaCDM + OmegaB;
  double Omega_sum1 = -Omega_M*exp_m1 - 2*OmegaR*exp_m2 + 2*OmegaLambda*exp_p2;
  double Omega_sum2 = Omega_M*exp_m1 + OmegaR*exp_m2 + OmegaLambda*exp_p2;
  double Omega_sum3 = Omega_M*exp_m1 + 4*OmegaR*exp_m2 + 4*OmegaLambda*exp_p2;

  return get_H0()*( Omega_sum3/(2*sqrt(Omega_sum2)) - (Omega_sum1*Omega_sum1)/(4*sqrt(Omega_sum2*Omega_sum2*Omega_sum2)) );
}


double BackgroundCosmology::get_OmegaB(double x) const{ 
  double Htemp = get_H0()/H_of_x(x);
  return OmegaB*exp(-3*x)*Htemp*Htemp;
}


double BackgroundCosmology::get_OmegaR(double x) const{ 
  double Htemp = get_H0()/H_of_x(x);
  return OmegaR*exp(-4*x)*Htemp*Htemp;
}


double BackgroundCosmology::get_OmegaNu(double x) const{ 
  return 0.0;
}


double BackgroundCosmology::get_OmegaCDM(double x) const{ 
  double Htemp = get_H0()/H_of_x(x);
  return OmegaCDM*exp(-3*x)*Htemp*Htemp;
}


double BackgroundCosmology::get_OmegaLambda(double x) const{ 
  double Htemp = get_H0()/H_of_x(x);
  return OmegaLambda*Htemp*Htemp;
}

double BackgroundCosmology::get_OmegaK(double x) const{ 
  return 0.0;
}

double BackgroundCosmology::eta_of_x(double x) const{
  return eta_of_x_spline(x);
}

double BackgroundCosmology::get_H0() const{ 
  return H0; 
}


double BackgroundCosmology::get_h() const{ 
  return h; 
}

double BackgroundCosmology::get_Neff() const{ 
  return Neff; 
}

double BackgroundCosmology::get_TCMB() const{ 
  return TCMB; 
}

double BackgroundCosmology::get_rho_crit() const{
  return rho_c;
}

//====================================================
// Print out info about the class
//====================================================
void BackgroundCosmology::info() const{ 
  std::cout << "\n";
  std::cout << "Info about cosmology class:\n";
  std::cout << "OmegaB:      " << OmegaB      << "\n";
  std::cout << "OmegaCDM:    " << OmegaCDM    << "\n";
  std::cout << "OmegaLambda: " << OmegaLambda << "\n";
  std::cout << "OmegaK:      " << OmegaK      << "\n";
  std::cout << "OmegaNu:     " << OmegaNu     << "\n";
  std::cout << "OmegaR:      " << OmegaR      << "\n";
  std::cout << "Neff:        " << Neff        << "\n";
  std::cout << "h:           " << h           << "\n";
  std::cout << "TCMB:        " << TCMB        << "\n";
  std::cout << std::endl;
} 

//====================================================
// Output some data to file
//====================================================
void BackgroundCosmology::output(const std::string filename) const{
  const double x_min = -15.0;
  const double x_max =  4.0;
  const int    n_pts =  10000;
  
  Vector x_array = Utils::linspace(x_min, x_max, n_pts);

  std::ofstream fp(filename.c_str());
  auto print_data = [&] (const double x) {
    fp << x                  << " ";  // 0
    fp << eta_of_x(x)        << " ";  // 1
    fp << Hp_of_x(x)         << " ";  // 2
    fp << dHpdx_of_x(x)      << " ";  // 3
    fp << get_OmegaB(x)      << " ";  // 4
    fp << get_OmegaCDM(x)    << " ";  // 5
    fp << get_OmegaLambda(x) << " ";  // 6
    fp << get_OmegaR(x)      << " ";  // 7
    fp << get_OmegaNu(x)     << " ";  // 8
    fp << get_OmegaK(x)      << " ";  // 9
    fp << ddHpddx_of_x(x)    << " ";  // 10
    fp <<"\n";
  };
  std::for_each(x_array.begin(), x_array.end(), print_data);
}

