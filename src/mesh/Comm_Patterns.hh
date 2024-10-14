//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   mesh/Comm_Patterns.hh
 * \author H. Park
 * \date   Aug. 1 2023
 * \brief  Comm_Patterns header file.
 * \note   Copyright (C) 2023-2024 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#ifndef rtt_mesh_Comm_Patterns_hh
#define rtt_mesh_Comm_Patterns_hh

#include "Topology.hh"
#include <map>
#include <memory>
#include <utility>

// cSpell:ignore topo

namespace rtt_mesh {

//================================================================================================//
/*!
 * \class Comm_Patterns
 * \brief Determine boundary cell communication patterns.

 * The Comm_Patterns class determines the boundary cell relationships between processors in a
 * domain-decomposed (DD) topology.  See rtt_mesh::Topology for a definition of \b boundary \b
 * cells. The Comm_Patterns class gives the following information:

 * \arg what processors this processor shares boundary cells with;
 * \arg what local cells this processor contains that another processor needs information about
 * \arg the boundary cells that a remote processor will send data to

 * The Comm_Patterns class defines information that the processor needs for \e gather/scatter type
 * operations.  This type of information is only required in spatially domain decomposed
 * topologies. In these cases a DD_Topology is always defined. The information contained in this
 * class is closely related to the information contained in the DD_Topology class.  However, the
 * information contained here is related to very specific topology-like information. Also, parallel
 * communication is required to generate the information; thus it belongs in its own class.

 * Accessor functions to Comm_Patterns are in the form of constant iterators to std::map data.  As
 * such, the user can access the data by iterating through it.  The std::map::iterator is a
 * std::pair so the processor is accessed by iterator->first and the cell list by
 * iterator->second[i].  Because all accessor functions to this data in Comm_Patterns are const,
 * only map::const_iterator types are returned.  This ensures that clients cannot modify the
 * Comm_Patterns data through the iterators. All accessing of the Comm_Patterns must be done through
 * const_iterator; expressions of the type Comm_Patterns::iterator i = object.get_recv_begin() will
 * be flagged as an error at compile time.  However, the iterators themselves can be modified (ie
 * Comm_Patterns::const_iterator i; i++; is valid).  See the examples for usage.
 *
 * Example usage of Comm_Patterns class. The full DD mesh is a nine cell mesh (3 x 3 with structured
 * ordering) with cells 1 and 2 on processor 0; cells 3 and 4 on processor 1; cells 5 and 6 on
 * processor 2; and cells 7, 8, and 9 on processor 4.
 *
 * \image html src/mesh/autodoc/DD_Comm_Patterns.png "2 proc/8 cell decomposed test mesh."
 * \image html src/mesh/autodoc/DD_Mesh.png "4 proc/9 cell decomposed test mesh"
 */
//================================================================================================//

class Comm_Patterns {
public:
  // Typedefs.
  using sf_int = std::vector<int>;
  using proc_map = std::map<int, sf_int>;
  using proc_pair = std::pair<int, sf_int>;
  using iterator = proc_map::iterator;
  using const_iterator = proc_map::const_iterator;
  using SP_Topology = std::shared_ptr<Topology>;

private:
  // Map of processor/local cells that this processor will send out to remote processors to fill
  // their boundary cell fields.
  proc_map send_query_map;

  // Map of processor/boundary cells that a remote processor will fill on a receive.
  proc_map recv_query_map;

private:
  // Implementation
  void calc_recv_map_DD(SP_Topology const &topo);
  void calc_send_map_DD(SP_Topology const &topo);

  // Empty the query maps.
  void empty_query_maps();

public:
  //! Constructor.
  explicit Comm_Patterns() noexcept
      : send_query_map(proc_map()), recv_query_map(proc_map()) { /*empty*/
  }

  // Calculate communication patterns.
  void calc_patterns(SP_Topology const &topo);

  // Query to see of the patterns are set.
  explicit operator bool() const;

  // ACCESSORS

  //! Get number of processors from which data is received.
  NODISCARD size_t get_num_recv_procs() const { return recv_query_map.size(); }

  //! Get number of processors to which data is sent.
  NODISCARD size_t get_num_send_procs() const { return send_query_map.size(); }

  //! Get iterator to beginning of receive data map.
  NODISCARD const_iterator get_recv_begin() const { return recv_query_map.begin(); }

  //! Get iterator to end of receive data map.
  NODISCARD const_iterator get_recv_end() const { return recv_query_map.end(); }

  //! Get iterator to beginning of send data map.
  NODISCARD const_iterator get_send_begin() const { return send_query_map.begin(); }

  //! Get iterator to end of send data map.
  NODISCARD const_iterator get_send_end() const { return send_query_map.end(); }
};

} // end namespace rtt_mesh

#endif // rtt_mesh_Comm_Patterns_hh

//------------------------------------------------------------------------------------------------//
// end of mesh/Comm_Patterns.hh
//------------------------------------------------------------------------------------------------//
