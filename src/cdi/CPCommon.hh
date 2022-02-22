//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   cdi/CPCommon.hh
 * \author Kendra Long
 * \date   Mon Jan 19 13:41:01 2001
 * \brief  Datatypes needed in CPEloss classes
 * \note   Copyright (C) 2019-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#ifndef rtt_cdi_CPCommon_hh
#define rtt_cdi_CPCommon_hh

#include <cstdint>

namespace rtt_cdi {

namespace constants {
//! Number of models contained in rtt_cdi::CPModel.
unsigned int const num_CPModels(3);

//! "Effectively infinite" stopping power (keV sh^-1)
constexpr double max_eloss = 1.e15;
} // namespace constants

enum class CPModelAngleCutoff {
  NONE = 0,  /*!< no angle cutoff */
  TNBURN = 1 /*!< cutoff angle formula from "TN Burn Project: Minimum Impact parameter Scattering
                  Angle Cutoff/Nuclear Interactions */
};

enum class CPModelType {
  UNASSIGNED_ETYPE = 0, /*!< unassigned type; used as a placeholder before deciding type */
  ANALYTIC_ETYPE = 1,   /*!< an Analytic eloss model */
  TABULAR_ETYPE = 2     /*!< a Tabular eloss model */
};

// Basic data required to describe a charged particle
class CParticle {
private:
  int32_t zaid;
  double mass;

public:
  int32_t get_zaid() const { return zaid; }
  int32_t get_z() const {
    if (zaid != -1) {
      // Integer division to extract the ZZZ from ZAID format ZZZAAA
      return zaid / 1000;
    } else {
      // Electrons have special ZAID = -1
      return -1;
    }
  }
  int32_t get_a() const { return zaid % 1000; }
  double get_mass() const { return mass; }

  CParticle(const int32_t zaid_in, const double mass_in) : zaid(zaid_in), mass(mass_in) {}
};

} // namespace rtt_cdi

#endif

//------------------------------------------------------------------------------------------------//
// End cdi/CPCommon.hh
//------------------------------------------------------------------------------------------------//
