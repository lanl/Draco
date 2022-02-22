//--------------------------------------------*-C++-*---------------------------------------------//
/*! \file   PhysicalConstexprs.hh
 *  \author Ben R. Ryan <brryan@lanl.gov>
 *  \brief  Provide a single place where mathematical and physical constants are defined in an
 *          entirely constexpr fashion.
 *  \date   2019 Nov 25
 *  \note   Copyright (C) 2019-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#ifndef rtt_units_PhysicalConstexprs_hh
#define rtt_units_PhysicalConstexprs_hh

#include "MathConstants.hh"

namespace rtt_units {
struct SI {
  static constexpr double Length = 1.;      // meter
  static constexpr double Mass = 1.;        // kilogram
  static constexpr double Time = 1.;        // second
  static constexpr double Temperature = 1.; // Kelvin
  static constexpr double Current = 1.;     // Amp
  static constexpr double Charge = 1.;      // Coulomb
  static constexpr double Capacitance = 1.; // Farad
  static constexpr double Angle = 1.;       // Radian
  static constexpr double Quantity = 1.;    // Mole
};

struct CGS {
  static constexpr double Length = 1.e2;                    // centimeter
  static constexpr double Mass = 1.e3;                      // gram
  static constexpr double Time = 1.;                        // second
  static constexpr double Temperature = 1.;                 // Kelvin
  static constexpr double Current = 1.e-1;                  // Biot
  static constexpr double Charge = 2.997924580e9;           // Statcoulomb
  static constexpr double Capacitance = 8.9831483395497e11; // Statfarad
  static constexpr double Angle = 1.;                       // Radian
  static constexpr double Quantity = 1.;                    // Mole
};

struct CGSH {
  static constexpr double Length = 1.e2;                    // centimeter
  static constexpr double Mass = 1.e3;                      // gram
  static constexpr double Time = 1.e8;                      // shake
  static constexpr double Temperature = 8.61733238496e-8;   // keV
  static constexpr double Current = 1.e-1;                  // Biot
  static constexpr double Charge = 2.997924580e9;           // Statcoulomb
  static constexpr double Capacitance = 8.9831483395497e11; // Statfarad
  static constexpr double Angle = 1.;                       // Radian
  static constexpr double Quantity = 1.;                    // Mole
};

//================================================================================================//
/*!
 * \class PhysicalConstexprs
 *
 * \brief A class to define and encapsulate physical and mathematical constants in a purely
 *        constexpr way.
 *
 * \example test/tstPhysicalConstexprs.cc
 * This is the unit regression test for the PhysicalConstexprs class.  It demonstrates typical
 * usage.
 */
//================================================================================================//
template <typename UNITS> class PhysicalConstexprs {
public:
  //! accesses Avogadro's number
  constexpr double avogadro() const { return d_avogadro; }
  //! see avogadro()
  constexpr double Na() const { return avogadro(); }

  //! accesses the fine structure constant
  constexpr double fineStructure() const { return d_fineStructure; }
  //! see fineStructure()
  constexpr double alpha() const { return fineStructure(); }

  //! accesses pi
  constexpr double pi() const { return d_pi; }

  //! accesses Planck's constant
  constexpr double planck() const { return d_planck; }
  //! see planck()
  constexpr double h() const { return planck(); }

  //! accesses the molar gas constant
  constexpr double gasConstant() const { return d_gasConstant; }
  //! see gasConstant()
  constexpr double R() const { return gasConstant(); }

  //! accesses the Boltzmann constant
  constexpr double boltzmann() const { return d_boltzmann; }
  //! see boltzmann()
  constexpr double k() const { return boltzmann(); }

  //! accesses the electron charge
  constexpr double electronCharge() const { return d_electronCharge; }
  //! see electronCharge()
  constexpr double e() const { return electronCharge(); }

  //! accesses the speed of light
  constexpr double speedOfLight() const { return d_speed_of_light; }
  //! see speedOfLight()
  constexpr double c() const { return d_speed_of_light; }

  //! accesses the Stefan-Boltzmann constant
  constexpr double stefanBoltzmann() const { return d_stefanBoltzmann; }
  //! see stefanBoltzmann()
  constexpr double sigma() const { return stefanBoltzmann(); }

  //! accesses the gravitational constant
  constexpr double gravitationalConstant() const { return d_gravitationalConstant; }
  //! see gravitationalConstant()
  constexpr double G() const { return gravitationalConstant(); }

  //! accesses the standard acceleration of gravity
  constexpr double accelerationFromGravity() const { return d_accelerationFromGravity; }
  //! see accelerationFromGravity()
  constexpr double g() const { return accelerationFromGravity(); }

  //! accesses the Faraday constant
  constexpr double faradayConstant() const { return d_faradayConstant; }
  //! see faradayConstant()
  constexpr double F() const { return faradayConstant(); }

  //! access the permeability of free space
  constexpr double permeabilityOfVacuum() const { return d_permeabilityOfVacuum; }
  //! see permeabilityOfVacuum()
  constexpr double mu0() const { return permeabilityOfVacuum(); }

  //! access the permittivity of free space
  constexpr double permittivityOfVacuum() const { return d_permittivityOfVacuum; }
  //! see permittivityOfVacuum()
  constexpr double epsi0() const { return permittivityOfVacuum(); }

  //! access the classical electron radius
  constexpr double classicalElectronRadius() const { return d_classicalElectronRadius; }
  //! see classicalElectronRadius()
  constexpr double re() const { return classicalElectronRadius(); }

  //! access the electron mass
  constexpr double electronMass() const { return d_electronMass; }
  //! see electronMass()
  constexpr double me() const { return electronMass(); }

  //! access the proton mass
  constexpr double protonMass() const { return d_protonMass; }
  //! see protonMass()
  constexpr double mp() const { return protonMass(); }

  //! access the Electronvolt
  constexpr double electronVolt() const { return d_electronVolt; }
  //! see electronVolt()
  constexpr double eV() const { return electronVolt(); }

  //! access the atomic mass unit
  constexpr double atomicMassUnit() const { return d_atomicMassUnit; }
  //! see atomicMassUnit()
  constexpr double amu() const { return atomicMassUnit(); }

private:
  //! Derived unit conversions
  static constexpr double Length = UNITS::Length;
  static constexpr double Mass = UNITS::Mass;
  static constexpr double Time = UNITS::Time;
  static constexpr double Temperature = UNITS::Temperature;
  static constexpr double Current = UNITS::Current;
  static constexpr double Charge = UNITS::Charge;
  static constexpr double Capacitance = UNITS::Capacitance;
  static constexpr double Angle = UNITS::Angle;
  static constexpr double Quantity = UNITS::Quantity;
  static constexpr double Force = Mass * Length / (Time * Time);
  static constexpr double Energy = Force * Length;
  static constexpr double Power = Energy / Time;

  //! Dimensionless constants
  //! Avogadro cosntant (CODATA 2010 value)
  static constexpr double d_avogadro = 6.02214129e23;

  //! Fine structure constant (CODATA 2010 value)
  static constexpr double d_fineStructure = 7.2973525698e-3;

  //! Pi
  static constexpr double d_pi = PI;

  //! SI values of fundamental constants
  //! Planck constant (CODATA 2010 value)
  static constexpr double d_planck_SI = 6.62606957e-34; // J s

  //! Molar gas constant (CODATA 2010 value)
  static constexpr double d_gasConstant_SI = 8.3144621; // J mol^-1 K^-1

  //! Boltzmann constant (CODATA 2010 value)
  static constexpr double d_boltzmann_SI = 1.380648800e-23; // J K^-1

  //! Electron charge (CODATA 2018 exact value)
  static constexpr double d_electronCharge_SI = 1.602176565e-19; // C

  //! Speed of light (CODATA 2018 exact value)
  static constexpr double d_speed_of_light_SI = 2.99792458e8; // m s^-1

  //! Gravitational constant (CODATA 2010 value)
  static constexpr double d_gravitationalConstant_SI = 6.67384e-11; // m^3 kg^-1 s^-2

  //! Standard acceleration of gravity (CODATA 2010 value)
  static constexpr double d_accelerationFromGravity_SI = 9.80665; // m s^-2

  //! Electron rest mass (CODATA 2010 value)
  static constexpr double d_electronMass_SI = 9.10938291e-31; // kg

  //! Proton rest mass (CODATA 2010 value)
  static constexpr double d_protonMass_SI = 1.672621777e-27; // kg

  //! Stefan-Boltzmann constant
  static constexpr double d_stefanBoltzmann_SI =
      2.0 * PI * PI * PI * PI * PI * d_boltzmann_SI * d_boltzmann_SI * d_boltzmann_SI *
      d_boltzmann_SI /
      (15.0 * d_planck_SI * d_planck_SI * d_planck_SI * d_speed_of_light_SI * d_speed_of_light_SI);

  //! Faraday constant
  static constexpr double d_faradayConstant_SI = d_avogadro * d_electronCharge_SI; // C mol^-1

  //! Permeability of free space
  static constexpr double d_permeabilityOfVacuum_SI = 4.0 * PI * 1.0e-7; // N A^-2

  //! Permittivity of free space
  static constexpr double d_permittivityOfVacuum_SI =
      1.0 / d_permeabilityOfVacuum_SI / d_speed_of_light_SI / d_speed_of_light_SI; // F m^-1

  //! Classical electron radius
  static constexpr double d_classicalElectronRadius_SI =
      d_electronCharge_SI * d_electronCharge_SI /
      (4 * PI * d_permittivityOfVacuum_SI * d_electronMass_SI * d_speed_of_light_SI *
       d_speed_of_light_SI); // m

  //! Electron volt
  static constexpr double d_electronVolt_SI = d_electronCharge_SI; // J

  //! Atomic mass unit (CODATA 2010 value)
  static constexpr double d_atomicMassUnit_SI = 1.660538921e-27; // kg

  //! Dimensional values converted to templated unit system
  //! Planck constant value in unit system
  static constexpr double d_planck = d_planck_SI * Energy * Time;

  //! Molar gas constant value in unit system
  static constexpr double d_gasConstant = d_gasConstant_SI * Energy / Quantity / Temperature;

  //! Boltzmann constant value in unit system
  static constexpr double d_boltzmann = d_boltzmann_SI * Energy / Temperature;

  //! Electron charge value in unit system
  static constexpr double d_electronCharge = d_electronCharge_SI * Charge;

  //! Speed of light value in unit system
  static constexpr double d_speed_of_light = d_speed_of_light_SI * Length / Time;

  //! Stefan-Boltzmann constant value in unit system
  static constexpr double d_stefanBoltzmann =
      d_stefanBoltzmann_SI * Power /
      (Length * Length * Temperature * Temperature * Temperature * Temperature);

  //! Gravitational constant value in unit system
  static constexpr double d_gravitationalConstant =
      d_gravitationalConstant_SI * Length * Length * Length / (Mass * Time * Time);

  //! Standard acceleration of gravity value in unit system
  static constexpr double d_accelerationFromGravity =
      d_accelerationFromGravity_SI * Length / (Time * Time);

  //! Faraday constant value in unit system
  static constexpr double d_faradayConstant = d_faradayConstant_SI * Capacitance;

  //! Permeability of vacuum value in unit system
  static constexpr double d_permeabilityOfVacuum =
      d_permeabilityOfVacuum_SI * Force / (Current * Current);

  //! Permittivity of vacuum value in unit system
  static constexpr double d_permittivityOfVacuum = d_permittivityOfVacuum_SI * Capacitance / Length;

  //! Classical electron radius value in unit system
  static constexpr double d_classicalElectronRadius = d_classicalElectronRadius_SI * Length;

  //! Electron mass value in unit system
  static constexpr double d_electronMass = d_electronMass_SI * Mass;

  //! Proton mass value in unit system
  static constexpr double d_protonMass = d_protonMass_SI * Mass;

  //! Electron volt value in unit system
  static constexpr double d_electronVolt = d_electronVolt_SI * Energy;

  //! Atomic mass unit value in unit system
  static constexpr double d_atomicMassUnit = d_atomicMassUnit_SI * Mass;
};
} // namespace rtt_units

#endif // rtt_units_PhysicalConstexprs_hh

//------------------------------------------------------------------------------------------------//
// End units/FundamentalConstexprs.hh
//------------------------------------------------------------------------------------------------//
