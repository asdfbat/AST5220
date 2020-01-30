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
  H0_SI = H0*1000/Constants.Mpc;
  OmegaR = pow(M_PI, 3)/15.0*pow(Constants.k_b*TCMB, 4)/(pow(Constants.hbar, 3)*pow(Constants.c, 5))*8*M_PI*Constants.G/(3*H0_SI*H0_SI);
  OmegaNu = 0;
  OmegaK = 0;
  printf("Sum of Omegas = %e\n", OmegaR + OmegaCDM + OmegaB + OmegaLambda);
  printf("OmegaR0 = %e\n", OmegaR);
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
  double x_stop = log(a_stop);
  printf("Solving conformal time for x in [%f, %f]\n", x_start, x_stop);

  Vector x_array = Utils::linspace(x_start, x_end, npts);


  // The ODE for deta/dx
  ODEFunction detadx = [&](double x, const double *eta, double *detadx){

    //=============================================================================
    // TODO: Set the rhs of the detadx ODE
    //=============================================================================
    //...
    //...

    detadx[0] = Constants.c/Hp_of_x(x);

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

  Vector eta_inc = {0.0};

  ODESolver ode;
  ode.solve(detadx, x_array, eta_inc);
  Vector2D all_data = ode.get_data();
  Vector eta_data(all_data.size());
  for(int i=0; i<all_data.size(); i++){
    eta_data[i] = all_data[i][0];
  }

  eta_of_x_spline.create(x_array, eta_data, "eta");

  Utils::EndTiming("Eta");
}


//====================================================
// Get methods
//====================================================


double BackgroundCosmology::H_of_x(double x) const{
  double Omega_sum = OmegaCDM*exp(-3*x) + OmegaB*exp(-3*x) + OmegaR*exp(-4*x) + OmegaLambda;
  // double Omega_sum = get_OmegaCDM(x) + get_OmegaB(x) + get_OmegaLambda(x) + get_OmegaR(x);
  return get_H0()*sqrt(Omega_sum);
}


double BackgroundCosmology::Hp_of_x(double x) const{
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
  double Htemp = H0/H_of_x(x);
  return OmegaB*exp(-3*x)*Htemp*Htemp;
}


double BackgroundCosmology::get_OmegaR(double x) const{ 
  double Htemp = H0/H_of_x(x);
  return OmegaR*exp(-4*x)*Htemp*Htemp;
}


double BackgroundCosmology::get_OmegaNu(double x) const{ 
  return 0.0;
}


double BackgroundCosmology::get_OmegaCDM(double x) const{ 
  double Htemp = H0/H_of_x(x);
  return OmegaCDM*exp(-3*x)*Htemp*Htemp;
}


double BackgroundCosmology::get_OmegaLambda(double x) const{ 
  double Htemp = H0/H_of_x(x);
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
    fp <<"\n";
  };
  std::for_each(x_array.begin(), x_array.end(), print_data);
}

