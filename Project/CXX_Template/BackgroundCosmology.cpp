#include "BackgroundCosmology.h"

//====================================================
// Constructors
//====================================================
    
BackgroundCosmology::BackgroundCosmology(
    double h, 
    double OmegaB, 
    double OmegaCDM, 
    double OmegaLambda,
    double Neff, 
    double TCMB) :
  h(h),
  OmegaB(OmegaB),
  OmegaCDM(OmegaCDM),
  OmegaLambda(OmegaLambda),
  Neff(Neff), 
  TCMB(TCMB)
{

  //=============================================================================
  // TODO: Compute OmegaR, OmegaNu, OmegaK, H0, ...
  //=============================================================================
  //...
  //...
  //...
  //...
  H0 = 100*h;
  OmegaR = pow(M_PI, 3)/15.0*pow(Constants.k_b*TCMB, 4)/(pow(Constants.hbar, 3)*pow(Constants.c, 4))*8*M_PI*Constants.G/(3*H0);
  OmegaNu = 0;
  OmegaK = 0;
  printf("Sum of Omegas = %e\n", OmegaR + OmegaCDM + OmegaB + OmegaLambda);
}

//====================================================
// Do all the solving. Compute eta(x)
//====================================================

// Solve the background
void BackgroundCosmology::solve(){
  Utils::StartTiming("Eta");
    
  //=============================================================================
  // TODO: Set the range of x and the number of points for the splines
  // For this Utils::linspace(x_start, x_end, npts) is useful
  //=============================================================================
  int npts = 1000;
  double a_start = 1e-7;
  double a_stop = 1;
  double x_start = log(a_start);
  double x_stop = log(a_start);

  Vector x_array = Utils::linspace(x_start, x_end, npts);


  // The ODE for deta/dx
  ODEFunction detadx = [&](double x, const double *eta, double *detadx){

    //=============================================================================
    // TODO: Set the rhs of the detadx ODE
    //=============================================================================
    //...
    //...

    detadx[0] = 0.0;

    return GSL_SUCCESS;
  };

  //=============================================================================
  // TODO: Set the initial condition, set up the ODE system, solve and make
  // the spline eta_of_x_spline 
  //=============================================================================
  // ...
  // ...
  // ...
  // ...

  Utils::EndTiming("Eta");
}

//====================================================
// Get methods
//====================================================

double BackgroundCosmology::H_of_x(double x) const{

  //=============================================================================
  // TODO: Implement...
  //=============================================================================
  //...
  //...

  double Omega_sum = get_OmegaCDM(x) + get_OmegaB(x) + get_OmegaLambda(x) + get_OmegaR(x);
  return get_H0()*sqrt(Omega_sum);
}

double BackgroundCosmology::Hp_of_x(double x) const{

  //=============================================================================
  // TODO: Implement...
  //=============================================================================
  //...
  //...

  return exp(x)*H_of_x(x);
}

double BackgroundCosmology::dHpdx_of_x(double x) const{

  //=============================================================================
  // TODO: Implement...
  //=============================================================================
  //...
  //...

  return 0.0;
}

double BackgroundCosmology::ddHpddx_of_x(double x) const{

  //=============================================================================
  // TODO: Implement...
  //=============================================================================
  //...
  //...

  return 0.0;
}

double BackgroundCosmology::get_OmegaB(double x) const{ 
  if(x == 0.0) return OmegaB;

  //=============================================================================
  // TODO: Implement...
  //=============================================================================
  //...
  //...

  return OmegaB*exp(-3*x);
}

double BackgroundCosmology::get_OmegaR(double x) const{ 
  if(x == 0.0) return OmegaR;

  //=============================================================================
  // TODO: Implement...
  //=============================================================================
  //...
  //...

  return OmegaR*exp(-4*x);
}

double BackgroundCosmology::get_OmegaNu(double x) const{ 
  if(x == 0.0) return OmegaNu;

  //=============================================================================
  // TODO: Implement...
  //=============================================================================
  //...
  //...

  return 0.0;
}

double BackgroundCosmology::get_OmegaCDM(double x) const{ 
  if(x == 0.0) return OmegaCDM;

  //=============================================================================
  // TODO: Implement...
  //=============================================================================
  //...
  //...

  return OmegaCDM*exp(-3*x);
}

double BackgroundCosmology::get_OmegaLambda(double x) const{ 
  if(x == 0.0) return OmegaLambda;

  //=============================================================================
  // TODO: Implement...
  //=============================================================================
  //...
  //...

  return OmegaLambda;
}

double BackgroundCosmology::get_OmegaK(double x) const{ 
  if(x == 0.0) return OmegaK;

  //=============================================================================
  // TODO: Implement...
  //=============================================================================
  //...
  //...

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
  const double x_min = -10.0;
  const double x_max =  0.0;
  const int    n_pts =  100;
  
  Vector x_array = Utils::linspace(x_min, x_max, n_pts);

  std::ofstream fp(filename.c_str());
  auto print_data = [&] (const double x) {
    fp << x                  << " ";
    fp << eta_of_x(x)        << " ";
    fp << Hp_of_x(x)         << " ";
    fp << dHpdx_of_x(x)      << " ";
    fp << get_OmegaB(x)      << " ";
    fp << get_OmegaCDM(x)    << " ";
    fp << get_OmegaLambda(x) << " ";
    fp << get_OmegaR(x)      << " ";
    fp << get_OmegaNu(x)     << " ";
    fp << get_OmegaK(x)      << " ";
    fp <<"\n";
  };
  std::for_each(x_array.begin(), x_array.end(), print_data);
}

