//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   RTT_Format_Reader/test/TestRTTPolyhedron.cc
 * \author James S. Warsa
 * \date   Fri Aug 19 12:31:47 2016
 * \brief  Testing the POLYHEDRON Element_Definition
 * \note   Copyright (C) 2016-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "RTT_Format_Reader/RTT_Mesh_Reader.hh"
#include "ds++/Release.hh"
#include "ds++/ScalarUnitTest.hh"
#include <sstream>

using rtt_mesh_element::Element_Definition;
using rtt_RTT_Format_Reader::RTT_Mesh_Reader;

using namespace std;
using namespace rtt_dsxx;
using namespace rtt_RTT_Format_Reader;

class is_cell {
public:
  //! Construct the predicate for a specified dimensionality.
  explicit is_cell(size_t const dimensionality_in) : dimensionality(dimensionality_in) { /* empty */
  }

  /*! \brief Returns \c true if the specified type index has the predicate dimensionality; \c false
   *         otherwise.
   */
  bool operator()(Element_Definition::Element_Type const type) {
    return (Element_Definition(type).get_dimension() == dimensionality);
  }

private:
  size_t dimensionality;
};

//------------------------------------------------------------------------------------------------//
// TESTS
//------------------------------------------------------------------------------------------------//
void test_polyhedron(rtt_dsxx::UnitTest &ut) {

  // Find the mesh file
  string const inpPath = ut.getTestSourcePath();
  vector<string> filenames = {"rttquad5.mesh", "rttquad9.mesh", "rttquad.mesh"};

  for (auto const &filename : filenames) {
    string const fq_filename(inpPath + filename);
    shared_ptr<RTT_Mesh_Reader> mesh(new RTT_Mesh_Reader(fq_filename));

    ostringstream m;
    m << "Read mesh file " << fq_filename << std::endl;
    ut.passes(m.str());

    size_t const ndim = mesh->get_dims_ndim();
    if (ndim != 2)
      FAILMSG("Unexpected dimension.");

    vector<shared_ptr<Element_Definition>> const element_defs(mesh->get_element_defs());
    for (size_t j = 0; j < element_defs.size(); ++j) {
      cout << "Element definition for element " << j << endl;
      element_defs[j]->print(cout);
    }
  }

  //----------------------------------------------------------------------------------------------//
  // Read polyhedron mesh file - this is the mesh that is of most interest in this test
  {
    vector<string> more_filenames = {"rttpolyhedron.mesh", "rttpolyhedron.2.mesh",
                                     "rttpolyhedron.2o.mesh", "rttpolyhedron.3.mesh",
                                     "rttpolyhedron.4.mesh"};

    for (unsigned i = 0; i < more_filenames.size(); ++i) {
      string filename(inpPath + more_filenames[i]);
      shared_ptr<RTT_Mesh_Reader> mesh(new RTT_Mesh_Reader(filename));

      ostringstream m;
      m << "Read mesh file " << filename << std::endl;
      ut.passes(m.str());

      // Investigate and report on the mesh

      // The element types begins with side types, followed by cell types.  We can distinguish these
      // by their dimensionality. Cell types have the full dimensionality of the mesh; we assume
      // side types have one less than the full dimensionality of the mesh.

      size_t const ndim = mesh->get_dims_ndim();
      vector<Element_Definition::Element_Type> const element_types(mesh->get_element_types());
      vector<vector<unsigned>> const element_nodes(mesh->get_element_nodes());
      map<string, set<unsigned>> element_sets(mesh->get_element_sets());
      vector<shared_ptr<Element_Definition>> const element_defs(mesh->get_element_defs());

      if (ndim != 3) {
        FAILMSG("Unexpected dimension.");
      } else {
        PASSMSG("Correct dimension.");
      }

      for (size_t j = 0; j < element_defs.size(); ++j) {
        cout << "Element definition for element " << j << endl;
        element_defs[j]->print(cout);
      }

      // typeof(auto) =
      // iterator_traits<vector<Element_Definition::Element_Type>>::difference_type
      auto const ncells = count_if(element_types.begin(), element_types.end(), is_cell(ndim));
      size_t const mcells = mesh->get_dims_ncells();
      Check(ncells >= 0);
      Check(mcells < INT_MAX);
      if (static_cast<size_t>(ncells) != mcells) {
        FAILMSG("Unexpected number of cells.");
      } else {
        PASSMSG("Correct number of cells.");
        std::cout << " There are " << ncells << " cells in the mesh" << std::endl;
      }

      Check(element_types.size() - ncells < UINT_MAX);
      auto nsides = static_cast<unsigned>(element_types.size() - ncells);
      size_t const msides = mesh->get_dims_nsides();
      if (nsides != msides) {
        FAILMSG("Unexpected number of sides.");
      } else {
        PASSMSG("Correct number of sides.");
        std::cout << " There are " << nsides << " sides in the mesh" << std::endl;
      }

      for (map<string, set<unsigned>>::const_iterator it = element_sets.begin();
           it != element_sets.end(); ++it) {
        // See if this is a side or cell flag.

        unsigned const representative_element = *it->second.begin();

        auto const type_index =
            static_cast<Element_Definition::Element_Type>(element_types[representative_element]);

        shared_ptr<Element_Definition const> type(new Element_Definition(type_index));

        if (type->get_dimension() == ndim) {
          std::cout << " Elements with flags " + it->first << " are cell elements " << std::endl;
        } else {
          // It has to be a side flag.
          if (type->get_dimension() + 1 != ndim) {
            throw invalid_argument(" Elements with flags " + it->first +
                                   " have the wrong dimension for a side element");
          } else {
            std::cout << " Elements with flags " + it->first << " are side elements" << std::endl;
          }
        }
      }

      unsigned ncorner = 0;
      for (unsigned et = nsides; et < element_types.size(); et++) {
        unsigned const c = et - nsides;
        shared_ptr<Element_Definition> cell_def = element_defs[c];
        unsigned const number_of_nodes = cell_def->get_number_of_nodes();
        ncorner += number_of_nodes;
      }

      if (i == 0 && ncorner != 13) {
        FAILMSG("Number of corners should be 13 for mesh 1");
      } else if (i == 1 && ncorner != 17) {
        FAILMSG("Number of corners should be 17 for mesh 2");
      } else {
        PASSMSG("Number of corners is correct");
      }

    } // End polyhedron mesh tests
  }
}

//------------------------------------------------------------------------------------------------//
int main(int argc, char *argv[]) {
  ScalarUnitTest ut(argc, argv, release);
  try {
    test_polyhedron(ut);
  }
  UT_EPILOG(ut);
}

//------------------------------------------------------------------------------------------------//
// end of TestRTTPolyhedron.cc
//------------------------------------------------------------------------------------------------//
