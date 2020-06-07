#include "Utils.h"
#include "BackgroundCosmology.h"
#include "RecombinationHistory.h"
#include "Perturbations.h"
#include "PowerSpectrum.h"

int main(int argc, char **argv){
  Utils::StartTiming("Everything");
  
  //=========================================================================
  // Parameters
  //=========================================================================

  // Background parameters
  double h           = 0.7;
  double OmegaB      = 0.046; // 0.05; //0.046;
  double OmegaCDM    = 0.224; //0.45;//0.224;
  double Neff        = 3.046;
  double TCMB        = 2.725;

  // Recombination parameters
  double Yp          = 0.0;

  //=========================================================================
  // Module I
  //=========================================================================
  // Set up and solve the background
  BackgroundCosmology cosmo(h, OmegaB, OmegaCDM, Neff, TCMB);
  cosmo.solve();
  cosmo.info();
  
  // Output background evolution quantities
  cosmo.output("../data/cosmology.txt");

  // Remove when module is completed
  // return 0;

  //=========================================================================
  // Module II
  //=========================================================================
  
  // Solve the recombination history
  RecombinationHistory rec(&cosmo, Yp);
  rec.solve();
  rec.info();
  rec.output("../data/recombination.txt");

  // Solve the recombination history using the Saha approximation the entire way.
  RecombinationHistory rec_saha(&cosmo, Yp);
  rec_saha.set_saha_limit(0.0);  // Setting the Xe limit of ending Saha regime to 0.
  rec_saha.solve();
  rec_saha.info();
  rec_saha.output("../data/recombination_saha.txt"); // Store data in a different file.

  // Remove when module is completed

  //=========================================================================
  // Module III
  //=========================================================================
 
  // Solve the perturbations
  Perturbations pert(&cosmo, &rec);
  pert.solve();
  pert.info();

  // Output perturbation quantities
  pert.output(1.0/Constants.Mpc, "../data/perturbations_k1.0.txt");
  pert.output(0.1/Constants.Mpc, "../data/perturbations_k0.1.txt");
  pert.output(0.01/Constants.Mpc, "../data/perturbations_k0.01.txt");
  pert.output(0.001/Constants.Mpc, "../data/perturbations_k0.001.txt");
  // return 0;
  
  //=========================================================================
  // Module IV
  //=========================================================================

  PowerSpectrum power(&cosmo, &rec, &pert);
  power.solve();
  power.output_Cell("../data/cells_Quad.txt");
  power.output_P("../data/P_Quad.txt");
  power.output_Theta("../data/Theta_l=2.txt", 0);
  power.output_Theta("../data/Theta_l=20.txt", 10);
  power.output_Theta("../data/Theta_l=120.txt", 20);
  power.output_Theta("../data/Theta_l=400.txt", 30);
  power.output_Theta("../data/Theta_l=900.txt", 40);
  power.output_Theta("../data/Theta_l=1400.txt", 50);
  power.output_Theta("../data/Theta_l=1900.txt", 60);
  
  // Remove when module is completed
  return 0;

  Utils::EndTiming("Everything");
}
