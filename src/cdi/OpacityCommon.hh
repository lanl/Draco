//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   cdi/OpacityCommon.hh
 * \author Kelly Thompson
 * \date   Mon Jan 19 13:41:01 2001
 * \brief  Datatypes needed in GrayOpacity and MultigroupOpacity
 * \note   Copyright (C) 2010-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#ifndef rtt_cdi_OpacityCommon_hh
#define rtt_cdi_OpacityCommon_hh

namespace rtt_cdi {

//================================================================================================//
// NUMBER OF MODELS AND REACTIONS
//================================================================================================//

namespace constants {
//! Number of models contained in rtt_cdi::Model.
unsigned int const num_Models(7);

//! Number of reaction types contained in rtt_cdi::Reaction.
unsigned int const num_Reactions(3);
} // namespace constants

//================================================================================================//
// ENUMERATIONS USED BY OPACITY CLASSES IN CDI
//================================================================================================//
/*!
 * \brief Physics model used to compute the opacity values.
 *
 * This enumeration \b must be unnumbered, ie it spans the set [0,N).  The number of models is given
 * by rtt_cdi::constants::num_Models.
 */
enum Model {
  ROSSELAND, /*!< use Rosseland mean opacities. */
  PLANCK,    /*!< use Plank mean opacities. */
  ANALYTIC,  /*!< use Analytic model opacities. */
  ISOTROPIC, /*!< use Isotropic scattering opacities. */
  THOMSON,   /*!< use Thomson scattering opacities. */
  COMPTON,   /*!< use Compton scattering opacities. */
  NOMODEL    /*!< null model */
};

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Opacity reaction type stored in this opacity object.
 *
 * This enumeration \b must be unnumbered, ie it spans the set [0,N).  The number of readtion types
 * is given by rtt_cdi::constants::num_Reactions.
 */
enum Reaction {
  TOTAL,      /*!< Total opacity value (scattering plus absorption). */
  ABSORPTION, /*!< Absorption cross sections only. */
  SCATTERING, /*!< Scattering cross sections only. */
  LAST_VALUE  /*!< dummy value */
};

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Type of opacity model: analytic, or gandolf.
 *
 * Used in Milagro Material_Data in packing the objects, returned by each opacity type. It was
 * previously defined as
 * \code
 *   typeid(rtt_cdi_analytic::Analytic_Odfmg_Opacity)
 * \endcode
 * mapping to 1, and
 * \code
 *   typeid(rtt_cdi_gandolf::GandolfOdfmgOpacity)
 * \endcode
 * mapping to 2.
 */
enum OpacityModelType {
  UNASSIGNED_TYPE = 0, /*!< unassigned type; used as a placeholder before deciding type */
  ANALYTIC_TYPE = 1,   /*!< an Analytic opacity model */
  GANDOLF_TYPE = 2,    /*!< a Gandolf opacity model */
  IPCRESS_TYPE = 3,    /*!< an Ipcress opacity model */
  DUMMY_TYPE = 99      /*!< a dummy opacity model */
};
} // end namespace rtt_cdi

#endif // rtt_cdi_OpacityCommon_hh

//------------------------------------------------------------------------------------------------//
// end of cdi/OpacityCommon.hh
//------------------------------------------------------------------------------------------------//
