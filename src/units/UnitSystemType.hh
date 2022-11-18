//--------------------------------------------*-C++-*---------------------------------------------//
/*! \file   UnitSystemType.hh
 *  \author Kelly Thompson
 *  \brief  Aggregates a collection of FundUnits to create a complete UnitSystemType.
 *  \date   Fri Oct 24 15:04:41 2003
 *  \note   Copyright (C) 2010-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#ifndef rtt_units_UnitSystemType_hh
#define rtt_units_UnitSystemType_hh

#include "FundUnit.hh"

namespace rtt_units {

//================================================================================================//
/*!
 * \class UnitSystemType
 *
 * \brief Aggregates a collection of FundUnits to create a complete UnitSystemType.  Also provides a
 *        tool for constructing UnitSystems on the fly.
 *
 * \sa UnitSystem
 * \sa test/tstUnitSystemType.cc
 * \sa test/tstUnitSystem.cc - a demonstration of using UnitSystemTypoe
 *
 * Different ways to construct a UnitSystem

 * \verbatim
 * using rtt_units::UnitSystemType;
 * using rtt_units::UnitSystem;
 * typedef rtt_units::UnitSystemType::
 *
 * UnitSystem uSI( UnitSystemType().SI() );
 * UnitSystem uX4( UnitSystemType().X4() );
 * UnitSystem uCGS( UnitSystemType().L( rtt_units::L_cm )
 *                                  .M( rtt_units::M_g )
 *                                  .t( rtt_units::t_s ) );
 * \endverbatim
 */
//================================================================================================//

class UnitSystemType {
public:
  // CONSTRUCTORS AND DESTRUCTOR

  //! default constructor
  UnitSystemType();
  ~UnitSystemType() = default;

  //! fully qualified constructor
  UnitSystemType(Ltype myL, Mtype myM, ttype myt, Ttype myT, Itype myI, Atype myA, Qtype myQ);

  //! Copy constructor
  UnitSystemType(UnitSystemType const &rhs);
  //! Move constructor
  UnitSystemType(UnitSystemType &&rhs) noexcept;

  //! Assignment
  UnitSystemType &operator=(UnitSystemType const &rhs) = default;
  //! Move Assignment
  UnitSystemType &operator=(UnitSystemType &&rhs) = default;

  // MANIPULATORS

  //! Set SI defaults
  UnitSystemType SI() { return {L_m, M_kg, t_s, T_K, I_amp, A_rad, Q_mol}; }

  //! Set X4 defaults
  UnitSystemType X4() { return {L_cm, M_g, t_sh, T_keV, I_amp, A_rad, Q_mol}; }

  //! Set cgs defaults
  UnitSystemType CGS() { return {L_cm, M_g, t_s, T_K, I_amp, A_rad, Q_mol}; }

  //! Set a FundUnit type for this UnitSystem

  UnitSystemType &L(Ltype myType, double const *cf = L_cf.data(),
                    std::string const &labels = L_labels);
  UnitSystemType &M(Mtype myType, double const *cf = M_cf.data(),
                    std::string const &labels = M_labels);
  UnitSystemType &t(ttype myType, double const *cf = t_cf.data(),
                    std::string const &labels = t_labels);
  UnitSystemType &T(Ttype myType, double const *cf = T_cf.data(),
                    std::string const &labels = T_labels);
  UnitSystemType &I(Itype myType, double const *cf = I_cf.data(),
                    std::string const &labels = I_labels);
  UnitSystemType &A(Atype myType, double const *cf = A_cf.data(),
                    std::string const &labels = A_labels);
  UnitSystemType &Q(Qtype myType, double const *cf = Q_cf.data(),
                    std::string const &labels = Q_labels);

  // ACCESSORS

  //! Return a FundUnit type when requested.

  FundUnit<Ltype> L() const { return d_L; }
  FundUnit<Mtype> M() const { return d_M; }
  FundUnit<ttype> t() const { return d_t; }
  FundUnit<Ttype> T() const { return d_T; }
  FundUnit<Itype> I() const { return d_I; }
  FundUnit<Atype> A() const { return d_A; }
  FundUnit<Qtype> Q() const { return d_Q; }

private:
  //! Fundamental unit types.

  FundUnit<Ltype> d_L;
  FundUnit<Mtype> d_M;
  FundUnit<ttype> d_t;
  FundUnit<Ttype> d_T;
  FundUnit<Itype> d_I;
  FundUnit<Atype> d_A;
  FundUnit<Qtype> d_Q;
};

} // end namespace rtt_units

#endif // rtt_units_UnitSystemType_hh

//------------------------------------------------------------------------------------------------//
// end of UnitSystemType.hh
//------------------------------------------------------------------------------------------------//
