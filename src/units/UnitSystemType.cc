//--------------------------------------------*-C++-*---------------------------------------------//
/*! \file   UnitSystemType.cc
 *  \author Kelly Thompson
 *  \brief  Aggregates a collection of FundUnits to create a complete UnitSystemType.
 *  \date   Fri Oct 24 15:04:41 2003
 *  \note   Copyright (C) 2010-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "UnitSystemType.hh"

namespace rtt_units {

//! default constructor
UnitSystemType::UnitSystemType()
    : d_L(FundUnit<Ltype>(L_null, L_cf.data(), L_labels)),
      d_M(FundUnit<Mtype>(M_null, M_cf.data(), M_labels)),
      d_t(FundUnit<ttype>(t_null, t_cf.data(), t_labels)),
      d_T(FundUnit<Ttype>(T_null, T_cf.data(), T_labels)),
      d_I(FundUnit<Itype>(I_null, I_cf.data(), I_labels)),
      d_A(FundUnit<Atype>(A_null, A_cf.data(), A_labels)),
      d_Q(FundUnit<Qtype>(Q_null, Q_cf.data(), Q_labels)) {}

//! fully qualified constructor
UnitSystemType::UnitSystemType(Ltype myL, Mtype myM, ttype myt, Ttype myT, Itype myI, Atype myA,
                               Qtype myQ)
    : d_L(FundUnit<Ltype>(myL, L_cf.data(), L_labels)),
      d_M(FundUnit<Mtype>(myM, M_cf.data(), M_labels)),
      d_t(FundUnit<ttype>(myt, t_cf.data(), t_labels)),
      d_T(FundUnit<Ttype>(myT, T_cf.data(), T_labels)),
      d_I(FundUnit<Itype>(myI, I_cf.data(), I_labels)),
      d_A(FundUnit<Atype>(myA, A_cf.data(), A_labels)),
      d_Q(FundUnit<Qtype>(myQ, Q_cf.data(), Q_labels)) {}

//! Copy constructor
UnitSystemType::UnitSystemType(UnitSystemType const &rhs)
    : d_L(rhs.L()), d_M(rhs.M()), d_t(rhs.t()), d_T(rhs.T()), d_I(rhs.I()), d_A(rhs.A()),
      d_Q(rhs.Q()) {}

//! Move constructor
UnitSystemType::UnitSystemType(UnitSystemType &&rhs) noexcept
    : d_L(rhs.L()), d_M(rhs.M()), d_t(rhs.t()), d_T(rhs.T()), d_I(rhs.I()), d_A(rhs.A()),
      d_Q(rhs.Q()) {}

//------------------------------------------------------------------------------------------------//

//! Set a FundUnit type for this UnitSystem

UnitSystemType &UnitSystemType::L(Ltype myType, double const *cf, std::string const &labels) {
  this->d_L = FundUnit<Ltype>(myType, cf, labels);
  return *this;
}
UnitSystemType &UnitSystemType::M(Mtype myType, double const *cf, std::string const &labels) {
  this->d_M = FundUnit<Mtype>(myType, cf, labels);
  return *this;
}
UnitSystemType &UnitSystemType::t(ttype myType, double const *cf, std::string const &labels) {
  this->d_t = FundUnit<ttype>(myType, cf, labels);
  return *this;
}
UnitSystemType &UnitSystemType::T(Ttype myType, double const *cf, std::string const &labels) {
  this->d_T = FundUnit<Ttype>(myType, cf, labels);
  return *this;
}
UnitSystemType &UnitSystemType::I(Itype myType, double const *cf, std::string const &labels) {
  this->d_I = FundUnit<Itype>(myType, cf, labels);
  return *this;
}
UnitSystemType &UnitSystemType::A(Atype myType, double const *cf, std::string const &labels) {
  this->d_A = FundUnit<Atype>(myType, cf, labels);
  return *this;
}
UnitSystemType &UnitSystemType::Q(Qtype myType, double const *cf, std::string const &labels) {
  this->d_Q = FundUnit<Qtype>(myType, cf, labels);
  return *this;
}

} // end namespace rtt_units

//------------------------------------------------------------------------------------------------//
// end of UnitSystemType.cc
//------------------------------------------------------------------------------------------------//
