//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   cdi_analytic/test/tstAnalytic_EICoupling.cc
 * \author Mathew Cleveland
 * \date   March 2019
 * \brief  Analytic_EICoupling test.
 * \note   Copyright (C) 2019-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "cdi/CDI.hh"
#include "cdi_analytic/Analytic_EICoupling.hh"
#include "ds++/Release.hh"
#include "ds++/ScalarUnitTest.hh"

using namespace std;

using rtt_cdi::CDI;
using rtt_cdi::EICoupling;
using rtt_cdi_analytic::Analytic_EICoupling;
using rtt_cdi_analytic::Analytic_EICoupling_Model;
using rtt_cdi_analytic::Constant_Analytic_EICoupling_Model;
using rtt_dsxx::soft_equiv;
using std::dynamic_pointer_cast;

//------------------------------------------------------------------------------------------------//
// TESTS
//------------------------------------------------------------------------------------------------//

void analytic_ei_coupling_test(rtt_dsxx::UnitTest &ut) {
  using Constant_Model = Constant_Analytic_EICoupling_Model;

  // make an analytic model (constant ei_coupling)
  shared_ptr<Constant_Model> model(new Constant_Model(1.1));

  if (!model)
    FAILMSG("Failed to build an Analytic EICoupling Model!");

  // make an analtyic electron-ion coupling
  Analytic_EICoupling analytic(model);

  // checks
  {
    double Te = 1.0;
    double Ti = 2.0;
    double rho = 3.0;

    double w_e = 4.0;
    double w_i = 5.0;
    double ei_coupling = 1.1;

    // specific heats
    if (!soft_equiv(analytic.getElectronIonCoupling(Te, Ti, rho, w_e, w_i), ei_coupling))
      ITFAILS;
  }

  // field check
  {
    vector<double> Te = {0.993, 0.882, 0.590, 0.112, 0.051, 0.001};
    vector<double> Ti = {0.993, 0.882, 0.590, 0.112, 0.051, 0.001};
    vector<double> rho = {0.993, 0.882, 0.590, 0.112, 0.051, 0.001};
    vector<double> w_e = {0.993, 0.882, 0.590, 0.112, 0.051, 0.001};
    vector<double> w_i = {0.993, 0.882, 0.590, 0.112, 0.051, 0.001};
    double ei_coupling_ref = 1.1;

    vector<double> ei_coupling = analytic.getElectronIonCoupling(Te, Ti, rho, w_e, w_i);

    if (ei_coupling.size() != 6)
      ITFAILS;

    for (int i = 0; i < 6; i++) {
      if (!soft_equiv(ei_coupling[i], ei_coupling_ref))
        ITFAILS;
    }
  }

  // Test the get_Analytic_Model() member function.
  {
    shared_ptr<Constant_Model const> myEICoupling_model =
        dynamic_pointer_cast<Constant_Model const>(analytic.get_Analytic_Model());
    shared_ptr<Constant_Model const> expected_model(model);

    if (expected_model == myEICoupling_model)
      PASSMSG("get_Analytic_Model() returned the expected EICoupling model.");
    else
      FAILMSG("get_Analytic_Model() did not return the expected EICoupling model.");
  }

  // Test the get_parameters() members function
  {
    std::vector<double> params(model->get_parameters());

    std::vector<double> expectedValue(1, 1.1);

    double const tol(1.0e-12);

    if (params.size() != expectedValue.size())
      ITFAILS;

    if (soft_equiv(params.begin(), params.end(), expectedValue.begin(), expectedValue.end(), tol))
      PASSMSG("get_parameters() returned the analytic expression coefficients.");
    else
      FAILMSG("get_parameters() did not return the analytic expression "
              "coefficients.");
  }

  return;
}

//------------------------------------------------------------------------------------------------//
void CDI_test(rtt_dsxx::UnitTest &ut) {
  using Constant_Model = Constant_Analytic_EICoupling_Model;

  // cdi object
  CDI eiCouplingData;

  // analytic model
  shared_ptr<Analytic_EICoupling_Model> model(new Constant_Model(1.1));

  // assign the electron-ion coupling object
  shared_ptr<Analytic_EICoupling> analytic_ei_coupling(new Analytic_EICoupling(model));

  // EICoupling object
  shared_ptr<const EICoupling> ei_coupling = analytic_ei_coupling;
  if (ei_coupling.get()) {
    auto &r = *ei_coupling.get();
    FAIL_IF_NOT(typeid(r) == typeid(Analytic_EICoupling));
  }

  // Assign the object to cdi
  eiCouplingData.setEICoupling(ei_coupling);

  // check
  if (!eiCouplingData.ei_coupling())
    FAILMSG("Can't reference EICoupling smart pointer");

  // make temperature and density fields
  vector<double> Te = {0.993, 0.882, 0.590, 0.112, 0.051, 0.001};
  vector<double> Ti = {0.993, 0.882, 0.590, 0.112, 0.051, 0.001};
  vector<double> rho = {0.993, 0.882, 0.590, 0.112, 0.051, 0.001};
  vector<double> w_e = {0.993, 0.882, 0.590, 0.112, 0.051, 0.001};
  vector<double> w_i = {0.993, 0.882, 0.590, 0.112, 0.051, 0.001};

  vector<double> EICoupling;

  // test the data
  {

    EICoupling = eiCouplingData.ei_coupling()->getElectronIonCoupling(Te, Ti, rho, w_e, w_i);

    if (EICoupling.size() != 6)
      ITFAILS;

    for (int i = 0; i < 6; i++) {
      double ei_coupling_ref = 1.1;
      if (!soft_equiv(EICoupling[i], ei_coupling_ref))
        ITFAILS;
    }
  }

  // now reset the CDI
  eiCouplingData.reset();

  // should catch this
  bool caught = false;
  try {
    eiCouplingData.ei_coupling();
  } catch (const rtt_dsxx::assertion & /* except */) {
    PASSMSG("Good, caught an unreferenced EICoupling shared_ptr!");
    caught = true;
  }
  if (!caught)
    FAILMSG("Failed to catch an unreferenced shared_ptr<EICoupling>!");

  // now assign the analytic electron-ion coupling to CDI directly
  eiCouplingData.setEICoupling(analytic_ei_coupling);
  FAIL_IF_NOT(eiCouplingData.ei_coupling());
  {
    auto sp = eiCouplingData.ei_coupling();
    if (sp.get()) {
      auto &r = *sp.get();
      FAIL_IF_NOT(typeid(r) == typeid(Analytic_EICoupling));
    } else {
      FAILMSG("unable to retrieve eiCouplingData.ei_coupling()/");
    }
  }

  // now test the data again

  // test the data
  {
    EICoupling = eiCouplingData.ei_coupling()->getElectronIonCoupling(Te, Ti, rho, w_e, w_i);

    if (EICoupling.size() != 6)
      ITFAILS;

    for (int i = 0; i < 6; i++) {
      double ei_coupling_ref = 1.1;
      if (!soft_equiv(EICoupling[i], ei_coupling_ref))
        ITFAILS;
    }
  }

  return;
}

//------------------------------------------------------------------------------------------------//
void packing_test(rtt_dsxx::UnitTest &ut) {
  using Constant_Model = Constant_Analytic_EICoupling_Model;

  vector<char> packed;
  {
    // make an analytic model (polynomial specific heats)
    shared_ptr<Constant_Model> model(new Constant_Model(1.1));

    // make an analtyic electron-ion coupling
    shared_ptr<EICoupling> ei_coupling(new Analytic_EICoupling(model));

    vector<char> internal_pack = ei_coupling->pack();
    packed.resize(internal_pack.size());
    packed = internal_pack;
  }

  Analytic_EICoupling n_ei_coupling(packed);

  // checks
  {
    double Te = 5.0;
    double Ti = 5.0;
    double rho = 3.0;
    double w_e = 3.0;
    double w_i = 3.0;

    double ei_coupling_ref = 1.1;

    // specific heats
    if (!soft_equiv(n_ei_coupling.getElectronIonCoupling(Te, Ti, rho, w_e, w_i), ei_coupling_ref))
      ITFAILS;
  }

  if (ut.numFails == 0)
    PASSMSG("Analytic EICoupling packing/unpacking test successfull.");

  return;
}

//------------------------------------------------------------------------------------------------//
int main(int argc, char *argv[]) {
  rtt_dsxx::ScalarUnitTest ut(argc, argv, rtt_dsxx::release);
  try {
    // >>> UNIT TESTS
    analytic_ei_coupling_test(ut);
    CDI_test(ut);
    packing_test(ut);
  }
  UT_EPILOG(ut);
}

//------------------------------------------------------------------------------------------------//
// end of tstAnalytic_EICoupling.cc
//------------------------------------------------------------------------------------------------//
