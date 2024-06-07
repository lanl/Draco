//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   special_functions/test/tstYlm.cc
 * \author Kent Budge
 * \date   Tue Jul  6 10:00:38 2004
 * \note   Copyright (C) 2010-2024 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "dsxx/Release.hh"
#include "dsxx/ScalarUnitTest.hh"
#include "dsxx/Soft_Equivalence.hh"
#include "special_functions/Ylm.hh"
#include "units/PhysicalConstants.hh"
#include <sstream>

using namespace std;
using namespace rtt_sf;

//------------------------------------------------------------------------------------------------//
// Tests
//------------------------------------------------------------------------------------------------//

void comparecPlk(unsigned const l, int const k, double const x, double const expVal,
                 rtt_dsxx::UnitTest &ut) {
  using rtt_dsxx::soft_equiv;
  std::ostringstream msg;

  double const val(cPlk(l, k, x));

  if (soft_equiv(val, expVal)) {
    msg << "cPlk(" << l << "," << k << ") function returned the expected value = " << expVal << ".";
    ut.passes(msg.str());
  } else {
    msg << "cPlk(" << l << "," << k << ") function did NOT return the expected value.\n"
        << "\tExpected " << expVal << ", but found " << val;
    ut.failure(msg.str());
  }
}

//------------------------------------------------------------------------------------------------//
void tstcPlk(rtt_dsxx::UnitTest &ut) {
  using rtt_dsxx::soft_equiv;
  using rtt_units::PI;

  cout << "\nTesting cPlk function.\n" << endl;

  std::ostringstream msg;
  double x(0.3568);

  // cPlk(0,0)
  double expVal(sqrt(1.0 / 4.0 / PI));
  comparecPlk(0, 0, x, expVal, ut);

  // cPlk(1,0)
  expVal = sqrt(3.0 / 4.0 / PI) * x;
  comparecPlk(1, 0, x, expVal, ut);

  // cPlk(1,1)
  expVal = -1.0 * sqrt(3.0 / 8.0 / PI) * sqrt(1.0 - x * x);
  comparecPlk(1, 1, x, expVal, ut);

  // cPlk(2,0)
  expVal = sqrt(5.0 / 4.0 / PI) / 2.0 * (3.0 * x * x - 1.0);
  comparecPlk(2, 0, x, expVal, ut);

  // cPlk(2,1)
  expVal = sqrt(5.0 / 24.0 / PI) * (-3.0 * x) * sqrt(1.0 - x * x);
  comparecPlk(2, 1, x, expVal, ut);

  // cPlk(2,2)
  expVal = sqrt(5.0 / 96.0 / PI) * 3.0 * (1.0 - x * x);
  comparecPlk(2, 2, x, expVal, ut);

  if (ut.dbcOn() && !ut.dbcNothrow()) {
    // Check out-of-bounds

    bool caught(false);
    try {
      comparecPlk(2, 3, x, expVal, ut);
    } catch (rtt_dsxx::assertion & /*err*/) {
      ut.passes("Caught out of bounds.");
      caught = true;
    }
    if (!caught)
      ut.failure("Did not catch out of bounds.");

    // Check mu out-of-range

    caught = false;
    x = -999999.999;
    try {
      comparecPlk(2, 0, x, expVal, ut);
    } catch (rtt_dsxx::assertion & /*err*/) {
      ut.passes("Caught mu out of range.");
      caught = true;
    }
    if (!caught)
      ut.failure("Did not catch mu out of range.");
  }

  return;
}

//------------------------------------------------------------------------------------------------//
void compareNormYlk(unsigned const l, int const k, double const theta, double const phi,
                    double const expVal, rtt_dsxx::UnitTest &ut) {
  using rtt_dsxx::soft_equiv;
  std::ostringstream msg;

  double const Ylk(normalizedYlk(l, k, theta, phi));

  if (soft_equiv(Ylk, expVal)) {
    msg << "normalizedYlk(" << l << "," << k
        << ") function returned the expected value = " << expVal << ".";
    ut.passes(msg.str());
  } else {
    msg << "normalizedYlk(" << l << "," << k << ") function did NOT return the expected value.\n"
        << "\tExpected " << expVal << ", but found " << Ylk;
    ut.failure(msg.str());
  }
}

//------------------------------------------------------------------------------------------------//
void compareRealYlk(unsigned const l, int const k, double const theta, double const phi,
                    double const expVal, rtt_dsxx::UnitTest &ut) {
  using rtt_dsxx::soft_equiv;
  std::ostringstream msg;

  double const Ylk(realYlk(l, k, theta, phi));

  if (soft_equiv(Ylk, expVal)) {
    msg << "realYlk(" << l << "," << k << ") function returned the expected value = " << expVal
        << ".";
    ut.passes(msg.str());
  } else {
    msg << "realYlk(" << l << "," << k << ") function did NOT return the expected value.\n"
        << "\tExpected " << expVal << ", but found " << Ylk;
    ut.failure(msg.str());
  }
}

//------------------------------------------------------------------------------------------------//
void compareGalerkinYlk(unsigned const l, int const k, double const theta, double const phi,
                        double const expVal, rtt_dsxx::UnitTest &ut) {
  using rtt_dsxx::soft_equiv;
  std::ostringstream msg;

  double const Ylk(galerkinYlk(l, k, cos(theta), phi, (4.0 * rtt_units::PI)));

  if (soft_equiv(Ylk, expVal)) {
    msg << "galerkinYlk(" << l << "," << k << ") function returned the expected value = " << expVal
        << ".";
    ut.passes(msg.str());
  } else {
    msg << "galerkinYlk(" << l << "," << k << ") function did NOT return the expected value.\n"
        << "\tExpected " << expVal << ", but found " << Ylk;
    ut.failure(msg.str());
  }
}

//------------------------------------------------------------------------------------------------//
void compareComplexYlk(unsigned const l, int const k, double const theta, double const phi,
                       double const expVal, rtt_dsxx::UnitTest &ut) {
  using rtt_dsxx::soft_equiv;
  std::ostringstream msg;

  double const Ylk(complexYlk(l, k, theta, phi));

  if (soft_equiv(Ylk, expVal)) {
    msg << "complexYlk(" << l << "," << k << ") function returned the expected value = " << expVal
        << ".";
    ut.passes(msg.str());
  } else {
    msg << "complexYlk(" << l << "," << k << ") function did NOT return the expected value.\n"
        << "\tExpected " << expVal << ", but found " << Ylk;
    ut.failure(msg.str());
  }
}

//------------------------------------------------------------------------------------------------//
void tstNormalizedYlk(rtt_dsxx::UnitTest &ut) {
  using rtt_dsxx::soft_equiv;
  using rtt_units::PI;

  cout << "\nTesting normalizedYlk function.\n" << endl;

  double const theta = 0.4, phi = 0.21;

  // y(0,0)
  double expVal = sqrt(1.0 / (4.0 * PI));
  compareNormYlk(0, 0, theta, phi, expVal, ut);

  // y(1,0)
  expVal = sqrt(3.0 / (4.0 * PI)) * cos(theta);
  compareNormYlk(1, 0, theta, phi, expVal, ut);
  // y(1,-1)
  expVal = sqrt(3.0 / (4.0 * PI)) * sin(theta) * cos(phi);
  compareNormYlk(1, -1, theta, phi, expVal, ut);
  // y(1,1)
  expVal = -1.0 * sqrt(3.0 / (4.0 * PI)) * sin(theta) * sin(phi);
  compareNormYlk(1, 1, theta, phi, expVal, ut);

  // y(2,-2)
  expVal = sqrt(15 / (16 * PI)) * sin(theta) * sin(theta) * cos(2 * phi);
  compareNormYlk(2, -2, theta, phi, expVal, ut);

  // y(2,-1)
  expVal = sqrt(15 / (4 * PI)) * sin(theta) * cos(theta) * cos(phi);
  compareNormYlk(2, -1, theta, phi, expVal, ut);

  // y(2,0)
  expVal = sqrt(5 / (16 * PI)) * (3 * cos(theta) * cos(theta) - 1);
  compareNormYlk(2, 0, theta, phi, expVal, ut);

  // y(2,1)
  expVal = -1.0 * sqrt(15 / (4 * PI)) * sin(theta) * cos(theta) * sin(phi);
  compareNormYlk(2, 1, theta, phi, expVal, ut);

  // y(2,2)
  expVal = sqrt(15 / (16 * PI)) * sin(theta) * sin(theta) * sin(2 * phi);
  compareNormYlk(2, 2, theta, phi, expVal, ut);

  // Test rotational invariance
  if (soft_equiv(normalizedYlk(3, 2, theta, phi), normalizedYlk(3, 2, theta + 2 * PI, phi)))
    ut.passes("Azimuthal invariance looks okay.");
  else
    ut.failure("Failed azimuthal invariance test.");
  if (soft_equiv(normalizedYlk(3, 2, theta, phi), normalizedYlk(3, 2, theta, phi + PI)))
    ut.passes("Polar invariance looks okay.");
  else
    ut.failure("Failed polar invariance test.");

  return;
}
//------------------------------------------------------------------------------------------------//
void tstRealYlk(rtt_dsxx::UnitTest &ut) {
  using rtt_dsxx::soft_equiv;
  using rtt_units::PI;

  cout << "\nTesting realYlk function.\n" << endl;

  double const theta = 0.9, phi = 0.06;

  // Y(0,0)
  double expVal = sqrt(1.0 / (4.0 * PI));
  compareRealYlk(0, 0, theta, phi, expVal, ut);

  // y(1,0)
  expVal = sqrt(3.0 / (4.0 * PI)) * cos(theta);
  compareRealYlk(1, 0, theta, phi, expVal, ut);
  // y(1,-1)
  double expVal_1_m1 = sqrt(3.0 / (8.0 * PI)) * sin(theta) * cos(phi);
  compareRealYlk(1, -1, theta, phi, expVal_1_m1, ut);
  // y(1,1)
  double expVal_1_1 = -1.0 * sqrt(3.0 / (8.0 * PI)) * sin(theta) * cos(phi);
  compareRealYlk(1, 1, theta, phi, expVal_1_1, ut);

  return;
}

//------------------------------------------------------------------------------------------------//
void tstComplexYlk(rtt_dsxx::UnitTest &ut) {
  using rtt_dsxx::soft_equiv;
  using rtt_units::PI;

  cout << "\nTesting complexYlk function.\n" << endl;

  double const theta = 2.96;
  double const phi = 4.27;

  // Y(0,0) only has a real component.
  double expVal = 0.0;
  compareComplexYlk(0, 0, theta, phi, expVal, ut);

  // y(1,0) only has a real component.
  expVal = 0.0;
  compareComplexYlk(1, 0, theta, phi, expVal, ut);
  // y(1,-1)
  double expVal_1_m1 = sqrt(3.0 / (8.0 * PI)) * sin(theta) * sin(phi);
  compareComplexYlk(1, -1, theta, phi, expVal_1_m1, ut);
  // y(1,1)
  double expVal_1_1 = -1.0 * sqrt(3.0 / (8.0 * PI)) * sin(theta) * sin(phi);
  compareComplexYlk(1, 1, theta, phi, expVal_1_1, ut);

  return;
}

//------------------------------------------------------------------------------------------------//
void tstgalerkinYlk(rtt_dsxx::UnitTest &ut) {
  using rtt_dsxx::soft_equiv;
  using rtt_units::PI;

  cout << "\nTesting realYlk function.\n" << endl;

  double const theta = 0.9, phi = 0.06;

  // Y(0,0)
  double expVal = (1.0 / (4.0 * PI));
  compareGalerkinYlk(0, 0, theta, phi, expVal, ut);

  //     // y(1,0)
  //     expVal = sqrt(3.0/(4.0*PI))*cos(theta);
  //     compareGalerkinYlk(1,0,theta,phi,expVal,ut);
  // y(1,-1)
  double expVal_1_m1 = -(3.0 / (4.0 * PI)) * sin(theta) * sin(phi);
  compareGalerkinYlk(1, -1, theta, phi, expVal_1_m1, ut);
  //     // y(1,1)
  double expVal_1_1 = -(3.0 / (4.0 * PI)) * sin(theta) * cos(phi);
  compareGalerkinYlk(1, 1, theta, phi, expVal_1_1, ut);

  return;
}

//------------------------------------------------------------------------------------------------//
int main(int argc, char *argv[]) {
  rtt_dsxx::ScalarUnitTest ut(argc, argv, rtt_dsxx::release);
  try {
    tstcPlk(ut);
    tstNormalizedYlk(ut);
    tstRealYlk(ut);
    tstComplexYlk(ut);
    tstgalerkinYlk(ut);
  }
  UT_EPILOG(ut);
}

//------------------------------------------------------------------------------------------------//
// end of testYlm.cc
//------------------------------------------------------------------------------------------------//
