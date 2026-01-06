//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   cdi_ndi/NDI_CP_Eloss.hh
 * \author Ben R. Ryan
 * \date   2020 Jun 3
 * \brief  NDI_CP_Eloss class definition.
 * \note   Copyright (C) 2020-2026 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#ifndef cdi_ndi_NDI_CP_Eloss_hh
#define cdi_ndi_NDI_CP_Eloss_hh

#define NDI_DEDX_SUPPORT (NDI_MAJOR > 2 || (NDI_MAJOR == 2 && NDI_MINOR >= 2))

#include "NDI_Base.hh"
#include "dsxx/Assert.hh"
#include "cdi/CPCommon.hh"
#include "cdi/CPEloss.hh"
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
  NDI_CP_Eloss(const std::string &library_in, const int zaidt_, const double mt_, const int zaidp_,
               const double mp_);

  //! Constructor (overridden gendir path)
  NDI_CP_Eloss(const std::string &gendir_in, const std::string &library_in, const int zaidt_,
               const double mt_, const int zaidp_, const double mp_);

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
  [[nodiscard]] double getEloss(const double temperature, const double density,
                                const double partSpeed,
                                const double debye_length = 0.0) const override;

  //! Query to see if data is in tabular or functional form (true)
  static constexpr bool is_data_in_tabular_form() { return true; }
  /*!
   * \brief Get target zaid
   */
  [[nodiscard]] int targ_zaid() const override { return zaidt; }

  /*!
   * \brief Get target mass
   */
  [[nodiscard]] double targ_mass() const override { return mt; }

  /*!
   * \brief Get projectile zaid
   */
  [[nodiscard]] int proj_zaid() const override { return zaidp; }

  /*!
   * \brief Get projectile mass
   */
  [[nodiscard]] double proj_mass() const override { return mp; }

  //! Get the material temperature grid as log(T) with units [keV]
  [[nodiscard]] sf_double getTemperatureGrid() const override { return temperatures; }

  //! Get the material density grid as log(p) with units [cm^{-3}]
  [[nodiscard]] sf_double getDensityGrid() const override { return densities; }

  //! Get the projectile energy grid as log(E) with units [MeV]
  [[nodiscard]] sf_double getEnergyGrid() const override { return energies; }

  //! Get the number of material temperature grid points.
  [[nodiscard]] size_t getNumTemperatures() const override { return n_temperature; }

  //! Get the number of material density grid points.
  [[nodiscard]] size_t getNumDensities() const override { return n_density; }

  //! Get the number of projectile energy grid points.
  [[nodiscard]] size_t getNumEnergies() const override { return n_energy; }

  //! Get the general eloss model type
  [[nodiscard]] rtt_cdi::CPModelType getModelType() const {
    return rtt_cdi::CPModelType::TABULAR_ETYPE;
  }

private:
// Only implemented if NDI is found
#ifdef NDI_FOUND
  void load_ndi();
#endif

private:
#ifdef NDI_FOUND
  rtt_units::PhysicalConstexprs<rtt_units::CGS> pc; //!< Units
#endif
  int zaidt{0};
  double mt{0.0};
  int zaidp{0};
  double mp{0.0};

  uint32_t n_energy{0};      //!< Number of gridpoints in projectile energy
  uint32_t n_density{0};     //!< Number of gridpoints in target density
  uint32_t n_temperature{0}; //!< Number of gridpoints in target temperature
#ifdef NDI_FOUND
  double d_log_energy{0.0};      //!< Log spacing of projectile energy gridpoints
  double d_log_density{0.0};     //!< Log spacing of target density gridpoints
  double d_log_temperature{0.0}; //!< Log spacing of target temperature gridpoints
  double min_energy{0.0};        //!< Minimum target energy
  double max_energy{0.0};        //!< Maximum target energy
  double min_density{0.0};       //!< Minimum target density
  double max_density{0.0};       //!< Maximum target density
  double min_temperature{0.0};   //!< Minimum target temperature
  double max_temperature{0.0};   //!< Maximum target temperature
#endif
  sf_double energies{};     //!< Vector of energy gridpoints
  sf_double densities{};    //!< Vector of density gridpoints
  sf_double temperatures{}; //!< Vector of temperature gridpoints
  // Note that after unit conversions, *_energy is really *_speed

  //! Storage for tabulated data
  sf_double stopping_data_1d{};

  //! 3D indexing of 1D stopping power data
  [[nodiscard]] inline double get_stopping_data(int i_t, int i_d, int i_e) const {
    return stopping_data_1d[i_e + n_energy * (i_d + n_density * i_t)];
  }
};

} // namespace rtt_cdi_ndi

#endif // cdi_ndi_NDI_CP_Eloss_hh

//------------------------------------------------------------------------------------------------//
// End cdi_ndi/NDI_CP_Eloss.hh
//------------------------------------------------------------------------------------------------//
