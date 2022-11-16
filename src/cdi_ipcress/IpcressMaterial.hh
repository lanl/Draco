//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   cdi_ipcress/IpcressMaterial.hh
 * \author Kelly Thompson
 * \date   Tue Aug 22 15:15:49 2000
 * \brief  Header file for IpcressMaterial class
 * \note   Copyright (C) 2011-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#ifndef rtt_cdi_ipcress_IpcressMaterial_hh
#define rtt_cdi_ipcress_IpcressMaterial_hh

#include "ds++/Assert.hh"
#include <algorithm>
#include <cctype>
#include <fstream>
#include <iostream>
#include <vector>

namespace rtt_cdi_ipcress {

//================================================================================================//
/*!
 * \class IpcressMaterial
 *
 * \brief Encapsulate the data associated with a single material loaded from an IPCRESS file.  This
 *        data includes:
 *
 * 1. The material identifier number (e.g.: 10001)
 *
 * 2. A list of data fields associated with the material.  This list will look something like this:
 *
 *    - tgrid
 *    - rgrid
 *    - hnugrid
 *    - ramg
 *    - rsmg
 *    - rgray
 *    - rtmg
 *    - pmg
 *    - pgray
 *    - tfree
 *    - tfree2
 *    - denoff
 *    - p
 *    - e
 *    - pelect
 *    - eelect
 *    - pnuc
 *    - enuc
 *    - hgnuplas
 *    - comp
 *    - frct
 *    - atwt
 *    - opname
 *    - general
 * .
 *
 * 3. The field values.  For example, 'tgrid' might be a list of 5 temperatures, {0.01, 0.1, 1.0,
 *    10.0, 100.0}.
 *
 * This object is used as a data container by the IpcressFile class.  A vector of empty
 * IpcressMaterials is created (length = num materials) and then the data for each material is
 * stored via the add_field() member function.
 */
//================================================================================================//

class IpcressMaterial {

  // NESTED CLASSES AND TYPEDEFS

  // DATA

  std::vector<std::string> fieldNames;
  std::vector<std::vector<double>> fieldValues;

  //! Ratio of Z to A for this material. Ration generated from 'comp' data.
  double zoa = {0.0};

public:
  // CREATORS

  //! Default constructor builds an empty object.
  IpcressMaterial() : fieldNames(), fieldValues() { /* empty */
  }

  // MANIPULATORS

  /*!
   * \brief Set the Z/A ratio for this material.
   * \param in_zoa Set the Z/A ratio for this material.
   */
  void set_zoa(double const in_zoa) { zoa = in_zoa; }

  /*!
   * \brief Add a field and it's data to the current material.  This is the normal way of populating
   *        this storage class. If a field is already registered, then overwrite the values.
   *
   * \param in_fieldName a string from the IPCRESS file that describes the associated data values
   *        (e.g.: tgrid, ramg, etc.)
   * \param in_values a vector<double> of values that represent the data loaded from the IPCRESS
   *        file.
   */
  void add_field(std::string &in_fieldName, std::vector<double> const &in_values) {
    // Remove white space from in_fieldName before saving it.
    // NOTE: ::isspace forces the use of c namespace rather than std::isspace
    in_fieldName.erase(std::remove_if(in_fieldName.begin(), in_fieldName.end(), ::isspace),
                       in_fieldName.end());

    // Check if the field name is already registered.
    auto itr = std::find(fieldNames.begin(), fieldNames.end(), in_fieldName);
    if (fieldNames.size() == 0 || itr == fieldNames.end()) {

      // Save the material data field (description and associated values).
      fieldNames.push_back(in_fieldName);
      fieldValues.push_back(in_values);

    } else {

      // Overwrite existing entries.
      size_t index = std::distance(fieldNames.begin(), itr);
      Check(fieldValues[index].size() == in_values.size());
      fieldValues[index] = in_values;
    }
    return;
  }

  // ACCESSORS

  //! return the list of known field descriptors
  std::vector<std::string> listDataFieldNames() const { return fieldNames; }

  //! return the vector of data associated with a field name.
  std::vector<double> data(std::string const &fieldName) const {
    Require(fieldName.size() > 0);
    Require(find(fieldNames.begin(), fieldNames.end(), fieldName) != fieldNames.end());
    return fieldValues[getFieldIndex(fieldName)];
  }

private:
  // IMPLEMENTATION

  //! Return the index of the provided string as stored in member data 'fieldNames'
  size_t getFieldIndex(std::string const &fieldName) const {
    Require(fieldName.size() > 0);
    Remember(auto pos = find(fieldNames.begin(), fieldNames.end(), fieldName));
    Check(pos != fieldNames.end());
    size_t fieldIndex = std::distance(fieldNames.begin(),
                                      std::find(fieldNames.begin(), fieldNames.end(), fieldName));
    Ensure(fieldIndex < fieldNames.size());
    return fieldIndex;
  }
};

} // end namespace rtt_cdi_ipcress

#endif // rtt_cdi_ipcress_IpcressFile_hh

//------------------------------------------------------------------------------------------------//
// end of cdi_ipcress/IpcressFile.hh
//------------------------------------------------------------------------------------------------//
