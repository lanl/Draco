//----------------------------------*-C++-*-----------------------------------//
/*!
 * \file   cdi_ndi/NDI_CP_Eloss.hh
 * \author Ben R. Ryan
 * \date   2020 Jun 3
 * \brief  NDI_CP_Eloss class definition.
 * \note   Copyright (C) 2019-2020 Triad National Security, LLC.
 *         All rights reserved. */
//----------------------------------------------------------------------------//

#ifndef cdi_ndi_NDI_CP_Eloss_hh
#define cdi_ndi_NDI_CP_Eloss_hh

#include "NDI_Base.hh"
#include "cdi/CPCommon.hh"
#include "cdi/CPEloss.hh"
#include "ds++/Assert.hh"
#include "units/PhysicalConstexprs.hh"

namespace rtt_cdi_ndi {

//============================================================================//
/*!
 * \class NDI_CP_Eloss
 *
 * \brief Derived rtt_cdi::CPEloss class for tabular eloss.
 *        This class implements the interface found in cdi/CPEloss.hh for
 *        the case where CP energy loss data is in tabular form, retrieved
 *        from NDI.
 */
//============================================================================//

class NDI_CP_Eloss : public rtt_cdi::CPEloss, public NDI_Base {
public:
  typedef std::vector<double> sf_double;

public:
  //! Constructor (default gendir path)
  NDI_CP_Eloss(const std::string &library_in, rtt_cdi::CParticle target_in,
               rtt_cdi::CParticle projectile_in);

  //! Constructor (overridden gendir path)
  NDI_CP_Eloss(const std::string &gendir_in, const std::string &library_in,
               rtt_cdi::CParticle target_in, rtt_cdi::CParticle projectile_in);

  //! Disable default constructor
  NDI_CP_Eloss() = delete;

  //! Disable copy constructor
  NDI_CP_Eloss(const NDI_CP_Eloss &) = delete;

  // >>> ACCESSORS

  //! Get a stopping power.
  double getEloss(const double temperature, const double density,
                  const double v0) const;

  //! Query to see if data is in tabular or functional form (true)
  static constexpr bool is_data_in_tabular_form() { return true; }

  //! Get the name of the associated data file
  inline std::string getDataFilename() const { return std::string(); }

  //! Get the material temperature grid.
  sf_double getTemperatureGrid() const { return temperatures; }

  //! Get the material density grid.
  sf_double getDensityGrid() const { return densities; }

  //! Get the projectile energy grid.
  sf_double getEnergyGrid() const { return energies; }

  //! Get the number of material temperature grid points.
  size_t getNumTemperatures() const { return n_temperature; }

  //! Get the number of material density grid points.
  size_t getNumDensities() const { return n_density; }

  //! Get the number of projectile energy grid points.
  size_t getNumEnergies() const { return n_energy; }

  //! Get the general eloss model type
  rtt_cdi::CPModelType getModelType() const {
    return rtt_cdi::CPModelType::TABULAR_ETYPE;
  }

private:
// Only implemented if NDI is found
#ifdef NDI_FOUND
  void load_ndi();
#endif

private:
  rtt_units::PhysicalConstexprs<rtt_units::CGS> pc; //!< Units

  uint32_t n_energy;          //!< Number of gridpoints in projectile energy
  uint32_t n_density;         //!< Number of gridpoints in target density
  uint32_t n_temperature;     //!< Number of gridpoints in target temperature
  double d_log_energy;        //!< Log spacing of projectile energy gridpoints
  double d_log_density;       //!< Log spacing of target density gridpoints
  double d_log_temperature;   //!< Log spacing of target temperature gridpoints
  double min_log_energy;      //!< Log of minimum projectile energy
  double min_log_density;     //!< Log of minimum target density
  double min_log_temperature; //!< Log of minimum target temperature
  double min_energy;          //!< Minimum target energy
  double max_energy;          //!< Maximum target energy
  double min_density;         //!< Minimum target density
  double max_density;         //!< Maximum target density
  double min_temperature;     //!< Minimum target temperature
  double max_temperature;     //!< Maximum target temperature
  sf_double energies;         //!< Vector of energy gridpoints
  sf_double densities;        //!< Vector of density gridpoints
  sf_double temperatures;     //!< Vector of temperature gridpoints
  // Note that after unit conversions, *_energy is really *_speed

  //! Storage for tabulated data
  sf_double stopping_data_1d;

  //! Path to gendir file, which indexes an NDI dataset
  std::string gendir;

  //! Name of library in which to find data
  const std::string library;
};

} // namespace rtt_cdi_ndi

#endif // cdi_ndi_NDI_CP_Eloss_hh

//----------------------------------------------------------------------------//
// End cdi_ndi/NDI_CP_Eloss.hh
//----------------------------------------------------------------------------//
