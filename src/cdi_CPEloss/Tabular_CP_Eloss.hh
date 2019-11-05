//----------------------------------*-C++-*----------------------------------//  
/*!                                                                              
 * \file   cdi_CPEloss/Tabular_CP_Eloss.hh                                     
 * \author Ben R. Ryan                                                        
 * \date   2019 Nov 4                                              
 * \brief  Tabular_CP_Eloss class definition.                                   
 * \note   Copyright (C) 2016-2019 Triad National Security, LLC.                 
 *         All rights reserved. */                                               
//---------------------------------------------------------------------------//

#ifndef __cdi_CPEloss_Tabular_CP_Eloss_hh__
#define __cdi_CPEloss_Tabular_CP_Eloss_hh__

#include "cdi/CPEloss.hh"
#include "ds++/DracoMath.hh"
#include <cmath>
#include <vector>

namespace rtt_cdi_cpeloss {
//===========================================================================//  
/*!                                                                              
 * \class Tabular_CP_Eloss                                                      
 *                                                                               
 * \brief Derived rtt_cdi::CPEloss class for tabular eloss.                     
 *        This class implements the interface found in cdi/CPEloss.hh for        
 *        the case where CP energy loss data is in tabular form, stored 
 *        in a file.                
 *                                                                               
 */                                                                              
//===========================================================================//

class Tabular_CP_Eloss : public rtt_cdi::CPEloss {
  public:
    // Useful typedefs
    using std::stod;
    using std::stoi;

  private:
  std::string filename;
  std::ifstream file;

  // Particle type being transported
  int32_t projectile_zaid;

  // Target species
  int32_t target_zaid;

  int32_t nbin_energy; // Number of bins in projectile energy
  int32_t nbin_density; // Number of bins in target density
  int32_t nbin_temperature; // Number of bins in target temperature
  double d_log_energy; // Width of projectile energy bin in log space
  double d_log_density; // Width of target density bin in log space
  double d_log_temperature; // Width of target temperature bin in log space
  double min_log_energy; // Log of minimum projectile energy
  double min_log_density; // Log of minimum target density
  double min_log_temperature; // Log of minimum target temperature

  // Utility for reading a line of an eloss file and as a vector of strings
  std::vector<std::string> read_line();

  public:
  // Constructor
  Tabular_CP_Eloss(std::string filename_in, int32_t target_zaid_in, int32_t projectile_zaid_in);

  double getEloss(const double temperature, const double density, const double v0) const;

  int32_t getTargetZAID() const { return target_zaid; }

  int32_t getProjectileZAID() const { return projectile_zaid; }
}

} // namespace rtt_cdi_cpeloss

#endif 