//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   cdi_ndi/NDI_CP_Eloss.hh
 * \author Ben R. Ryan
 * \date   2020 Jun 3
 * \brief  NDI_CP_Eloss class definition.
 * \note   Copyright (C) 2020-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#ifndef cdi_ndi_NDI_CP_Eloss_hh
#define cdi_ndi_NDI_CP_Eloss_hh

#define NDI_DEDX_SUPPORT (NDI_MAJOR > 2 || (NDI_MAJOR == 2 && NDI_MINOR >= 2))

#include "NDI_Base.hh"
#include "cdi/CPCommon.hh"
#include "cdi/CPEloss.hh"
#include "ds++/Assert.hh"
#include "units/PhysicalConstexprs.hh"

namespace rtt_cdi_ndi {

//================================================================================================//
/*!
 * \class NDI_CP_Eloss
 *
 * \brief Derived rtt_cdi::CPEloss class for tabular eloss.  This class implements the interface
 *        found in cdi/CPEloss.hh for the case where CP energy loss data is in tabular form,
 *        retrieved from NDI.
 */
//================================================================================================//

class NDI_CP_Eloss : public rtt_cdi::CPEloss, public NDI_Base {
public:
  using sf_double = std::vector<double>;

public:
  //! Constructor (default gendir path)
  NDI_CP_Eloss(const std::string &library_in, rtt_cdi::CParticle target_in,
               rtt_cdi::CParticle projectile_in);

  //! Constructor (overridden gendir path)
  NDI_CP_Eloss(const std::string &gendir_in, const std::string &library_in,
               rtt_cdi::CParticle target_in, rtt_cdi::CParticle projectile_in);

  //! Default destructor
  ~NDI_CP_Eloss() override = default;

  //! Disable default constructor
  NDI_CP_Eloss() = delete;

  //! Disable copy constructor
  NDI_CP_Eloss(const NDI_CP_Eloss &rhs) = delete;

  //! Disable move constructor
  NDI_CP_Eloss(NDI_CP_Eloss &&rhs) noexcept = delete;

  //! Disable assignment
  NDI_CP_Eloss &operator=(const NDI_CP_Eloss &rhs) = delete;

  //! Disable move assignment
  NDI_CP_Eloss &operator=(NDI_CP_Eloss &&rhs) noexcept = delete;

  // >>> ACCESSORS

  //! Get a stopping power.
  double getEloss(const double temperature, const double density,
                  const double partSpeed) const override;

  //! Query to see if data is in tabular or functional form (true)
  static constexpr bool is_data_in_tabular_form() { return true; }

  //! Get the material temperature grid.
  sf_double getTemperatureGrid() const override { return temperatures; }

  //! Get the material density grid.
  sf_double getDensityGrid() const override { return densities; }

  //! Get the projectile energy grid.
  sf_double getEnergyGrid() const override { return energies; }

  //! Get the number of material temperature grid points.
  size_t getNumTemperatures() const override { return n_temperature; }

  //! Get the number of material density grid points.
  size_t getNumDensities() const override { return n_density; }

  //! Get the number of projectile energy grid points.
  size_t getNumEnergies() const override { return n_energy; }

  //! Get the general eloss model type
  rtt_cdi::CPModelType getModelType() const { return rtt_cdi::CPModelType::TABULAR_ETYPE; }

private:
// Only implemented if NDI is found
#ifdef NDI_FOUND
  void load_ndi();
#endif

private:
#ifdef NDI_FOUND
  rtt_units::PhysicalConstexprs<rtt_units::CGS> pc; //!< Units
#endif

  uint32_t n_energy{0};      //!< Number of gridpoints in projectile energy
  uint32_t n_density{0};     //!< Number of gridpoints in target density
  uint32_t n_temperature{0}; //!< Number of gridpoints in target temperature
#ifdef NDI_FOUND
  double d_log_energy{0.0};        //!< Log spacing of projectile energy gridpoints
  double d_log_density{0.0};       //!< Log spacing of target density gridpoints
  double d_log_temperature{0.0};   //!< Log spacing of target temperature gridpoints
  double min_log_energy{0.0};      //!< Log of minimum projectile energy
  double min_log_density{0.0};     //!< Log of minimum target density
  double min_log_temperature{0.0}; //!< Log of minimum target temperature
  double min_energy{0.0};          //!< Minimum target energy
  double max_energy{0.0};          //!< Maximum target energy
  double min_density{0.0};         //!< Minimum target density
  double max_density{0.0};         //!< Maximum target density
  double min_temperature{0.0};     //!< Minimum target temperature
  double max_temperature{0.0};     //!< Maximum target temperature
#endif
  sf_double energies{};     //!< Vector of energy gridpoints
  sf_double densities{};    //!< Vector of density gridpoints
  sf_double temperatures{}; //!< Vector of temperature gridpoints
  // Note that after unit conversions, *_energy is really *_speed

  //! Storage for tabulated data
  sf_double stopping_data_1d{};

  //! 3D indexing of 1D stopping power data
  inline double get_stopping_data(int i_e, int i_d, int i_t) const {
    return stopping_data_1d[i_t + n_temperature * (i_d + n_density * i_e)];
  }
};

} // namespace rtt_cdi_ndi

#endif // cdi_ndi_NDI_CP_Eloss_hh

//------------------------------------------------------------------------------------------------//
// End cdi_ndi/NDI_CP_Eloss.hh
//------------------------------------------------------------------------------------------------//
