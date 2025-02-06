//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   mesh/Topology.hh
 * \author H. Park
 * \date   Aug. 1 2023
 * \brief  Topology class header file.
 * \note   Copyright (C) 2023-2025 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#ifndef rtt_mesh_Topology_hh
#define rtt_mesh_Topology_hh

#include "dsxx/config.h" // defines [[nodiscard]]
#include <iosfwd>
#include <string>
#include <vector>

namespace rtt_mesh {

//================================================================================================//
/*!
 * \class Topology
 *
 * \brief Topology is a parallel-topology abstract base class that gives the processor/cell layout.
 *
 * The Topology keeps a record of processor<-->mesh-cell relationships.  The Topology class contains
 * all the information necessary for a processor to communicate with other processors that share
 * common or adjacent parts of a mesh.  Consider the following 2 processor decomposed mesh:
 *
 * \image html DD_Comm_Patterns.png "2 processor decomposed mesh."
 *
 * The cells in blue are local to a processor while the cells in green are \b ghost \b cells.
 * Ghost cells are cells that a processor needs information about to complete its calculation.
 * The Topology base class
 * defines accessors for the following data:
 *
 * \param parallel_scheme a string that returns the type of parallel decomposition:
 *        - "replication" domain replication
 *        - "DD" domain decomposition
 * \param num_local_cells the number of cells on a processor.
 * \param num_cells the total number of cells in the problem.
 * \param global_cell the global cell index of a local cell.
 * \param local_cell the local cell index of a global cell that lives on the processor.
 * \param remote_local_cell the local cell index of a ghost cell defined on the processor where
 *          the ghost cell resides
 * \param remote_proc the processor where a boundary cell resides
 *
 * Topology has two derived classes:
 * - rtt_mesh::DD_Topology class that is used for domain decomposed topologies
 * - rtt_mesh::Rep_Topology class that is used for full replication topologies.
 *
 * The topology family of classes should be build by an appropriate builder class.  The constructors
 * for the derived classes demand different information.  See rtt_mesh::DD_Topology and
 * rtt_mesh::Rep_Topology for the appropriate building requirements.
 *
 * Local and Global cell indices run [1,N].  Zero is reserved for special purposes (when a queried
 * cell does not reside on processor).
 *
 * \sa rtt_mesh::DD_Topology and rtt_mesh::Rep_Topology
 */
//================================================================================================//
class Topology {
public:
  // >>> NESTED TYPES

  // Typedefs.
  using packed_type = std::vector<char>;

  /*!
   * \brief Parallel decomposition enumeration.
   *
   * This type enumeration is not implemented yet by most code; however, we will be going to this
   * enumeration as a replacement for the string lookups. */
  enum Decomposition {
    DOMAIN_REPLICATION = 0,   //!< Each PE owns all (overlapping)
    DOMAIN_DECOMPOSITION = 1, //!< Each PE owns a non-overlapping section
    NUM_DECOMPOSITION = 2
  };

protected:
  // >>> DATA

  //! Number of global cells.
  int64_t num_global_cells_d;

  //! Parallel scheme for this topology.
  std::string parallel_scheme_d;

public:
  // Constructor.
  explicit Topology(std::string scheme);

  // Constructor.
  Topology(const int64_t num_global_cells, std::string scheme);

  // copy constructor..
  Topology(const Topology &other_topo) = delete;
  // copy assignment
  Topology &operator=(const Topology &other_topo) = delete;
  // move ctor
  Topology(Topology &&) = delete;
  //move assignment
  Topology &operator=(Topology &&) = delete;

  //! Destructor for proper performance in inheritance tree.
  virtual ~Topology() = default;

  //! Return the parallel scheme.
  [[nodiscard]] const std::string &get_parallel_scheme() const { return parallel_scheme_d; }

  //! Return the number of global cells.
  [[nodiscard]] int64_t num_cells() const { return num_global_cells_d; }

  [[nodiscard]] bool valid_global_cell(int cell) const {
    return (cell >= 1) && (cell <= num_global_cells_d);
  }

  // >>> VIRTUAL FUNCTIONS SHARED BY DERIVED CLASSES

  //! Return the number of cells on processor.
  [[nodiscard]] virtual int64_t num_local_cells() const = 0;

  [[nodiscard]] virtual bool valid_local_cell(int cell) const = 0;

  //! Get the global cell index on-processor.
  [[nodiscard]] virtual int64_t global_cell(size_t local_cell) const = 0;

  //! Get the local cell index on-processor.
  [[nodiscard]] virtual int local_cell(int64_t global_cell) const = 0;

  //! Get the number of boundary cells on a processor.
  [[nodiscard]] virtual size_t num_ghost_cells() const = 0;

  //! Get the remote local cell index for a boundary cell.
  [[nodiscard]] virtual int remote_local_cell(size_t boundary_cell) const = 0;

  //! Get the remote processor id for a boundary cell.
  [[nodiscard]] virtual int remote_proc(size_t boundary_cell) const = 0;

  //! Get the global offset at a rank
  [[nodiscard]] virtual int64_t get_global_offset(int rank) const = 0;

  //! Packing function.
  [[nodiscard]] virtual packed_type pack() const = 0;

  //! Diagnostic printing.
  virtual void print(std::ostream &out) const = 0;
};

//------------------------------------------------------------------------------------------------//
// OVERLOADED OPERATORS
//------------------------------------------------------------------------------------------------//

//! Overloaded stream output operator.
std::ostream &operator<<(std::ostream &out, const Topology &topo);

} // end namespace rtt_mesh

#endif // rtt_mesh_Topology_hh

//------------------------------------------------------------------------------------------------//
// end of mesh/Topology.hh
//------------------------------------------------------------------------------------------------//
