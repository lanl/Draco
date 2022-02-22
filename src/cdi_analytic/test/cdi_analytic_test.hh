//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   cdi_analytic/test/cdi_analytic_test.hh
 * \author Thomas M. Evans
 * \date   Mon Sep 24 12:04:00 2001
 * \brief  Dummy model used for testing cdi_analytic software.
 * \note   Copyright (C) 2010-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#ifndef rtt_cdi_analytic_test_hh
#define rtt_cdi_analytic_test_hh

#include "cdi_analytic/Analytic_Models.hh"
#include "ds++/Packing_Utils.hh"

namespace rtt_cdi_analytic_test {

//================================================================================================//
// USER-DEFINED ANALYTIC_OPACITY_MODEL
//================================================================================================//

class Marshak_Model : public rtt_cdi_analytic::Analytic_Opacity_Model {
private:
  double a;

public:
  explicit Marshak_Model(double a_) noexcept : a(a_) { /*...*/
  }

  double calculate_opacity(double T, double /*rho*/) const override { return a / (T * T * T); }

  double calculate_opacity(double T, double rho, double /*nu*/) const override {
    return calculate_opacity(T, rho);
  }

  double calculate_opacity(double T, double rho, double /*nu0*/, double /*nu1*/) const override {
    return calculate_opacity(T, rho);
  }

  std::vector<double> get_parameters() const override { return {a}; }

  std::vector<char> pack() const override {
    rtt_dsxx::Packer packer;
    std::vector<char> p(sizeof(double) + sizeof(int));
    packer.set_buffer(p.size(), &p[0]);
    int indicator = 10;
    packer << indicator << a;
    return p;
  }
};

} // end namespace rtt_cdi_analytic_test

#endif // rtt_cdi_analytic_test_hh

//------------------------------------------------------------------------------------------------//
// end of cdi_analytic/test/cdi_analytic_test.hh
//------------------------------------------------------------------------------------------------//
