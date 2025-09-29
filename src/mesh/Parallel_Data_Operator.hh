//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   mesh/Parallel_Data_Operator.hh
 * \author H.Park
 * \date   Aug 1 2023
 * \brief  Parallel_Data_Operator header file.
 * \note   Copyright (C) 2023-2025 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#ifndef rtt_mesh_Parallel_Data_Operator_hh
#define rtt_mesh_Parallel_Data_Operator_hh

#include "Comm_Patterns.hh"
#include "c4/global.hh"

namespace rtt_mesh {

//================================================================================================//
/*!
 * \class Parallel_Data_Operator
 *
 * \brief Class to operate on field-type data across processors.
 *
 * Parallel_Data_Operator takes field-type arguments and performs operations on
 * the data across processor space.  The fields must have random access
 * iterators.  Thus, vectors, c-style arrays, deques, etc. will work.
 * Additional functionality is required for containers and fields that have
 * iterators with less functionality than random access.
 *
 * The Parallel_Data_Operator uses the rtt_mesh::Topology class to determine
 * mapping between processors.
 *
 * While most functions in this class are not inlined, they are contained in the
 * header file (template functions have internal linkage).  Thus, explicit
 * instantiations on different field types are not necessary.
 *
 * This class also contains two equivalence-type functions to check data across
 * processor space.  These are intended to be used for DBC checking.
 */
/*!
 *
 * Example usages of the Parallel_Data_Operator class.  In particular, note how
 * in functions test_mapping_DD() and test_mapping_replication() the nested
 * operation classes are used to determine the mapping for various data fields
 * in parallel topologies.
 *
 * \image html src/mesh/autodoc/DD_Mesh.png "4 proc/9 cell decomposed test mesh"
 */
//================================================================================================//

class Parallel_Data_Operator {
public:
  using SP_Topology = std::shared_ptr<Topology>;
  using SP_Comm_Patterns = std::shared_ptr<Comm_Patterns>;

  // Nested data operations classes.
  struct Data_Replicated;
  struct Data_Decomposed;
  struct Data_Distributed;

private:
  // Nested container classes for data for communication.
  template <typename T> class Data_Request;
  template <typename T> class Data_Send;

  //! Smart pointer to a topology class.
  SP_Topology topology;

public:
  //! Constructor.
  explicit Parallel_Data_Operator(SP_Topology top) : topology(top) { // NOLINT
    Ensure(topology != nullptr);
  }

  int64_t gid(size_t lid) { return topology->global_cell(lid); }
  [[nodiscard]] int get_remote_proc(size_t boundary_cell) const {
    return topology->remote_proc(boundary_cell);
  }
  [[nodiscard]] int get_remote_local_cell(size_t boundary_cell) const {
    return topology->remote_local_cell(boundary_cell);
  }
  [[nodiscard]] int64_t get_global_offset(int rank) const {
    return topology->get_global_offset(rank);
  }

  //! Map local data fields to global data fields.
  template <typename FT, typename GT, typename Op>
  void local_to_global(FT &local, GT &global, Op mapping);

  //! Gather boundary cell data for all processors
  template <typename FT, typename BT>
  void gather_ghost_cell_data(SP_Comm_Patterns const &pattern, const FT &local_data,
                              BT bc_data) const;

  //! Gather boundary cell data for all processors
  //! assume local_data and bc_data is vector.
  template <typename FT, typename BT>
  void gather_ghost_cell_data(SP_Comm_Patterns const &pattern, const FT &local_data, BT bc_data,
                              size_t const vector_size) const;

  //! Gather boundary cell data for all processors using vectors.
  template <typename T>
  void gather_ghost_cell_data(SP_Comm_Patterns pattern, const std::vector<T> &local_data,
                              std::vector<T> &bc_data);

  //! Gather boundary cell data for all processors using vectors.
  template <typename T>
  void gather_ghost_cell_data(SP_Comm_Patterns pattern,
                              const std::vector<std::vector<T>> &local_data,
                              std::vector<std::vector<T>> &bc_data, size_t const vector_size);

  //! Gather boundary cell data for all processors with per-ghost-cell variable vector lengths
  //! assume local_data and bc_data are vectors.
  template <typename T>
  std::vector<std::vector<T>>
  gather_ghost_cell_data(SP_Comm_Patterns const &pattern,
                         const std::vector<std::vector<T>> &local_data,
                         const std::vector<size_t> &local_vector_size,
                         const std::vector<size_t> &ghost_vector_size) const;

  //! Gather boundary cell data for all processors, where the data can be sparse (some cell data not
  //! populated), with per-ghost-cell variable vector lengths assume local_data and bc_data are
  //! vectors.
  template <typename T>
  std::vector<std::vector<T>>
  gather_sparse_ghost_cell_data(SP_Comm_Patterns const &pattern,
                                const std::vector<std::vector<T>> &local_data,
                                const std::vector<size_t> &local_vector_size,
                                const std::vector<size_t> &ghost_vector_size) const;

  //! Gather boundary cell vector data for all processors with per-ghost-cell variable lengths
  //! assume local_data and bc_data are vectors.
  template <typename T>
  std::vector<std::vector<std::vector<T>>>
  gather_ghost_cell_data(SP_Comm_Patterns const &pattern,
                         const std::vector<std::vector<std::vector<T>>> &local_data,
                         const std::vector<size_t> &local_vector_size,
                         const std::vector<size_t> &ghost_vector_size, const int Nv) const;

  //! Ghost cell data exchange
  template <typename T>
  void fill_ghost_cell_data(SP_Comm_Patterns pattern, std::vector<T> &local_data,
                            const size_t vector_size = 1);

protected:
  //! Gather boundary cell data for all processors
  //! assume local_data and bc_data is FLAT vector.
  template <typename FT, typename BT>
  void gather_fv_ghost_cell_data(SP_Comm_Patterns const &pattern, const FT &local_data, BT bc_data,
                                 size_t const vector_size) const;
};

} // namespace rtt_mesh

//------------------------------------------------------------------------------------------------//
// MEMBER TEMPLATE DEFINITIONS
//------------------------------------------------------------------------------------------------//

#include "Parallel_Data_Operator.i.hh"

#endif // rtt_mesh_Parallel_Data_Operator_hh

//------------------------------------------------------------------------------------------------//
// end of mesh/Parallel_Data_Operator.hh
//------------------------------------------------------------------------------------------------//
