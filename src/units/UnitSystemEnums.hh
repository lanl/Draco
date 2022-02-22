//--------------------------------------------*-C++-*---------------------------------------------//
/*! \file   UnitSystemEnums.hh
 *  \author Kelly Thompson
 *  \brief  This file contains enums, conversion factors and labels that help define a UnitSystem.
 *  \date   Fri Oct 24 15:57:09 2003
 *  \note   Copyright (C) 2010-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#ifndef rtt_units_UnitSystemEnums_hh
#define rtt_units_UnitSystemEnums_hh

#include "ds++/config.h"
#include <array>
#include <string>

namespace rtt_units {

//========================================//
// ENUMERATED LENGTH TYPES
//========================================//

enum Ltype {
  L_null = 0, //!< no length type
  L_m = 1,    //!< meters (m)
  L_cm = 2    //!< centimeters (1 m = 100 cm)
};

int constexpr num_Ltype = 3;
constexpr std::array<double, num_Ltype> L_cf = {0.0, 1.0, 100.0};
constexpr char const *L_labels = "NA,m,cm";
constexpr char const *L_long_labels = "no length unit specified,meter,centimeter";

//========================================//
// ENUMERATED MASS TYPES
//========================================//

enum Mtype {
  M_null = 0, //!< no mass type
  M_kg = 1,   //!< kilogram (kg)
  M_g = 2     //!< gram (1 g = 1e-3 kg)
};

int constexpr num_Mtype = 3;
constexpr std::array<double, num_Mtype> M_cf = {0.0, 1.0, 1000.0};
constexpr char const *M_labels = "NA,kg,g";
constexpr char const *M_long_labels = "no mass unit specified,kilogram,gram";

//========================================//
// ENUMERATED TIME TYPES
//========================================//

enum ttype {
  t_null, //!< no time type
  t_s,    //!< seconds (s)
  t_ms,   //!< milliseconds (1 ms = 1e-3 s)
  t_us,   //!< microseconds (1 us = 1e-6 s)
  t_sh,   //!< shakes       (1 ns = 1e-8 s)
  t_ns    //!< nanoseconds  (1 ns = 1e-9 s)
};

int constexpr num_ttype = 6;
constexpr std::array<double, num_ttype> t_cf = {0.0, 1.0, 1.0e3, 1.0e6, 1.0e8, 1.0e9};
constexpr char const *t_labels = "NA,s,ms,us,sh,ns";
constexpr char const *t_long_labels =
    "no time unit specified,second,milisecond,microsecond,shake,nanosecond";

//========================================//
// ENUMERATED TEMPERATURE TYPES
//========================================//

enum Ttype {
  T_null, //!< no temperature type
  T_K,    //!< Kelvin
  T_keV,  //!< keV (1 K=8.61733238496e-8 keV or 1 keV=1.1604519302808940e7 K). This conversion
          //!factor between K and keV must agree with the value given in PhysicalConstants.hh.
  T_eV    //!< eV (1 K = 8.61733238496e-5 keV or 1 eV = 11604.519302808940 K)
};

int constexpr num_Ttype = 4;
constexpr std::array<double, num_Ttype> T_cf = {0.0, 1.0, 1.0 / 1.1604519302808940e7,
                                                1.0 / 1.1604519302808940e4};
constexpr char const *T_labels = "NA,K,keV,eV";

//========================================//
// ENUMERATED CURRENT TYPES
//========================================//

enum Itype {
  I_null, //!< no current type0.001
  I_amp   //!< Amp (SI)
};

int constexpr num_Itype = 2;
constexpr std::array<double, num_Itype> I_cf = {0.0, 1.0};
constexpr char const *I_labels = "NA,Amp";

//========================================//
// ENUMERATED ANGLE TYPES
//========================================//

enum Atype {
  A_null, //!< no angle type
  A_rad,  //!< Radian (SI)
  A_deg   //!< degrees (PI rad = 180 deg)
};

int constexpr num_Atype = 3;
constexpr std::array<double, num_Atype> A_cf = {0.0, 1.0, 57.295779512896171};
constexpr char const *A_labels = "NA,rad,deg";

//========================================//
// ENUMERATED QUANTITY TYPES
//========================================//

enum Qtype {
  Q_null, //!< no quantity type
  Q_mol   //!< mole (SI)
};

int constexpr num_Qtype = 2;
constexpr std::array<double, num_Qtype> Q_cf = {0.0, 1.0};
constexpr char const *Q_labels = "NA,mol";

//------------------------------------------------------------------------------------------------//
// HELPER FUNCTIONS
//------------------------------------------------------------------------------------------------//

//! Extract unit labels from list in UnitSystemEnums.hh.
std::string setUnitLabel(size_t const pos, std::string const &labels);

} // end namespace rtt_units

#endif // rtt_units_UnitSystemEnums_hh

//------------------------------------------------------------------------------------------------//
// end of UnitSystemEnums.hh
//------------------------------------------------------------------------------------------------//
