//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   mesh/Parallel_Data_Operator.i.hh
 * \author H. Park
 * \date   Aug 1 2023
 * \brief  Member definitions of class Parallel_Data_Operator
 * \note   Copyright (C) 2023-2025 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#ifndef mesh_Parallel_Data_Operator_i_hh
#define mesh_Parallel_Data_Operator_i_hh

#include <iterator>
#include <limits>

// cSpell:ignore LCOV, gsum, itor

namespace rtt_mesh {

//------------------------------------------------------------------------------------------------//
// LOCAL-TO-GLOBAL MAPPING OPERATIONS
//------------------------------------------------------------------------------------------------//
/*!
 * \brief Map local data fields to global data fields.
 *
 * This function takes a local field and maps the data to a global field across all processors.  The
 * function is templated on field type (FT) and Op type.  The Op type provides a functor class
 * (Data_Decomposed, Data_Distributed, Data_Replicated) that describes how data should be mapped
 * across processors in certain topologies.  The following list describes the data combinations:
 *
 * \arg \b "replication"/Data_Replicated local_field is copied directly into the global_field, it is
 *      assumed that the local_fields are the same on each processor;
 *
 * \arg \b "replication"/Data_Decomposed the local_field is copied into the global field on each
 *      processor, the global_fields are subsequently summed across all processors;
 *
 * \arg \b "DD"/Data_Distributed local_fields are mapped into the appropriate cells on the
 *      global_field on each processor, the global_fields are then summed across all processors;
 *
 * \arg \b "DD/replication"/Data_Decomposed first the local_field data is mapped into local copies
 *      of the global_field on each processor, the global_fields are then summed-up across all
 *      processors;
 *
 * \arg \b "DD/replication"/Data_Replicated first the local_field data is mapped into local copies
 *      of the global_field on each processor, the global_fields are then mapped across all
 *      processors with replicated cells being assigned the proper value and non-replicated cells
 *      communicating its results to all other processors.
 *
 * The Parallel_Data_Operator knows the problem topology.  The user must specify an appropriate
 * nested operations functor class when calling this function.  The choices are
 * Parallel_Data_Operator::Data_Replicated, Parallel_Data_Operator::Data_Decomposed, and
 * Parallel_Data_Operator::Data_Distributed.  See the examples for info.
 *
 * \param local_field local field data on a processor
 * \param global_field mutable global field data sized to the number of global cells.
 * \param mapping nested operations functor class (see above)
 */
template <typename FT, typename GT, typename Op>
void Parallel_Data_Operator::local_to_global(FT &local_field, GT &global_field, Op mapping) {
  // the global field should be equal to the number of global cells and exists on every processor
  Require(static_cast<int>(global_field.size()) == static_cast<int>(topology->num_cells()));
  Require(global_field.size() > 0);

  // if topology is replication we do a simple copy
  if (topology->get_parallel_scheme() == "replication") {
    // the local field and global field should be the same size for replication topologies
    Check(local_field.size() == global_field.size());

    // perform the local-global mapping depending upon the data-type
    mapping(local_field, global_field, *this);

  } else if (topology->get_parallel_scheme() == "DD") {
    Check(static_cast<int>(local_field.size()) == static_cast<int>(topology->num_local_cells()));
    Check(static_cast<int>(global_field.size()) == static_cast<int>(topology->num_cells()));

    // iterators for local field
    auto begin = local_field.begin();
    auto end = local_field.end();

    // iterators for global field
    auto global_begin = global_field.begin();
    Remember(auto global_end = global_field.end());

    // iterators for reading/writing
    typename FT::iterator itr;
    typename GT::iterator global;

    // sweep through local cells and add them to the appropriate point in the local processor
    // version of the global field
    for (itr = begin; itr != end; itr++) {
      // calculate the local cell index
      auto local_cell = std::distance(begin, itr) + 1;
      Check(local_cell > 0 && static_cast<int>(local_cell) <= static_cast<int>(local_field.size()));

      // calculate the global cell index
      auto global_cell = topology->global_cell(local_cell);

      // calculate the global cell iterator position corresponding to the global_cell index -> the
      // advance function is much faster when global is a random access iterator
      global = global_begin;
      std::advance(global, (global_cell - 1));
      Check(global >= global_begin && global < global_end);

      // assign the local cell value to the global cell value
      *global = *itr;
    }

    // now do a summation over all processors of the local versions of the global fields
    Check(global_field.size() < INT32_MAX);
    rtt_c4::global_sum(&global_field[0], static_cast<int>(global_field.size()));
  } else {
    Insist(0, "Haven't done the rest yet!"); // LCOV_EXCL_LINE
  }
}

//------------------------------------------------------------------------------------------------//
// NESTED OPERATIONS CLASSES
//------------------------------------------------------------------------------------------------//

//! Nested operations functor class for data-replicated data.
struct Parallel_Data_Operator::Data_Replicated {
  // operator for performing data replicated operations in full replication topologies
  template <typename FT>
  void operator()(FT &local_field, FT &global_field, Parallel_Data_Operator const & /*data_op*/) {
    Require(local_field.size() == global_field.size());

    auto begin = local_field.begin();
    auto end = local_field.end();
    auto global = global_field.begin();

    // on each processor assign local data into global data
    std::copy(begin, end, global);
  }

  template <typename FT>
  void operator()(FT & /**/, FT & /**/, SP_Topology const & /**/,
                  Parallel_Data_Operator const & /**/) {
    // this will be a combination of data_distributed and data_replicated operations

    // first do a data distributed (full DD-like) operation locally

    // second cannot do a blind global sum: for cells that only live on one processor, do a gsum(),
    // for cells that are replicated, assign if they are on the processor, communicate to processors
    // that don't have the cell --> then check on the global_field
  }
};

//------------------------------------------------------------------------------------------------//
//! Nested operations functor class for data-decomposed data.
struct Parallel_Data_Operator::Data_Decomposed {
  // operator for performing data decomposed operations in full replication topologies
  template <typename FT>
  void operator()(FT &local_field, FT &global_field, Parallel_Data_Operator const & /*data_op*/) {
    Require(local_field.size() == global_field.size());

    auto local_begin = local_field.begin();
    auto local_end = local_field.end();
    auto global_begin = global_field.begin();

    // on each processor assign local data into global data
    std::copy(local_begin, local_end, global_begin);

    // sum up global data
    Check(global_field.size() < INT32_MAX);
    rtt_c4::global_sum(&global_field[0], static_cast<int>(global_field.size()));
  }

  template <typename FT, typename GT>
  void operator()(FT & /*unused*/, GT & /*unused*/, SP_Topology const & /*unused*/,
                  Parallel_Data_Operator const & /*unused*/) {
    // this will be a combination of data_distributed and data_decomposed operations

    // first do a data distributed (full DD-like) operation locally

    // second do a blind global_sum
  }
};

//------------------------------------------------------------------------------------------------//
//! Nested operations functor class for data-distributed data.
struct Parallel_Data_Operator::Data_Distributed {
  template <typename FT, typename GT>
  void operator()(FT & /*local*/, GT & /*global*/, Parallel_Data_Operator const & /*data_op*/) {
    // Do nothing
  }
};

//------------------------------------------------------------------------------------------------//
//! A Nested template class which contains memory for data and the C4 object for receiving it.
template <typename T> class Parallel_Data_Operator::Data_Request {
  T *data{nullptr};
  size_t data_size{0};
  rtt_c4::C4_Req request;

public:
  Data_Request() : request() { /* ... */
  }
  explicit Data_Request(size_t size) { allocate(size); }

  ~Data_Request() { release(); }

  // move ctor
  Data_Request(Data_Request &&) = delete;
  // move assignment
  Data_Request &operator=(Data_Request &&) = delete;

  void release() {
    delete[] data;
    data = nullptr;
    data_size = 0;
  }
  void allocate(size_t size) {
    T *temp = new T[size];
    assign(temp, size);
  }

  // C4_Req.complete should be const
  [[nodiscard]] bool complete() { return request.complete(); }
  [[nodiscard]] bool inuse() const { return request.inuse(); }

  void post(int processor, int tag) {
    Check(data_size < INT32_MAX);
    rtt_c4::receive_async(request, data, static_cast<int>(data_size), processor, tag);
  }

  T *get_data() { return data; }
  [[nodiscard]] T get_data(int c) const { return data[c]; }

  [[nodiscard]] size_t size() const { return data_size; }

  //! Copy constructor
  Data_Request(Data_Request const &rhs) : data(nullptr), data_size(rhs.data_size), request() {
    // This is used when constructing a vector of Data_Request objects.
    data_size = rhs.data_size;
    allocate(data_size);
    std::copy(rhs.data, rhs.data + data_size, data);
    return;
  }

  // Disable assignment operator
  Data_Request &operator=(Data_Request const & /*rhs*/);

private:
  void assign(T *data_, size_t size) {
    release();
    data = data_;
    data_size = size;
  }
};

//------------------------------------------------------------------------------------------------//
//! A Nested Template class which contains memory and implements the C4 operations for sending it
template <typename T> class Parallel_Data_Operator::Data_Send {
  T *data{nullptr};
  size_t data_size{0};

public:
  Data_Send() = default;
  explicit Data_Send(size_t size) { allocate(size); }
  ~Data_Send() { release(); }

  // Disable copy construction and assignment operators
  Data_Send(Data_Send const & /*rhs*/);
  Data_Send &operator=(Data_Send const & /*rhs*/);

  // move ctor
  Data_Send(Data_Send &&) = delete;
  // move assignment
  Data_Send &operator=(Data_Send &&) = delete;

  void release() {
    delete[] data;
    data = nullptr;
    data_size = 0;
  }
  void allocate(size_t size) {
    T *temp = new T[size];
    assign(temp, size);
  }

  void send(int processor, int tag) {
    Check(data_size < INT32_MAX);
    rtt_c4::send(data, static_cast<int>(data_size), processor, tag);
  }

  T *get_data() { return data; }
  T &get_data(int c) { return data[c]; }

  [[nodiscard]] size_t size() const { return data_size; }

private:
  void assign(T *data_, size_t size) {
    release();
    data = data_;
    data_size = size;
  }
};

//------------------------------------------------------------------------------------------------//
// GHOST CELL VECTOR PARALLEL OPERATIONS (GATHER/SCATTER TYPE)
//------------------------------------------------------------------------------------------------//
/*!
 * \brief Fill boundary cell vector data fields in a spatially decomposed topology.
 *
 * This function utilizes the rtt_mesh::Comm_Patterns class to "gather" data into ghost cell
 * fields.  What this entails is the following:
 *
 * \arg using the Comm_Patterns, each processor posts requests to processors it needs data from;

 * \arg using the Comm_Patterns and local data fields, each processor sends data out to processors
 * that need data;

 * \arg each processor receives the data for its boundary cells and writes it into a local boundary
 * cell field

 * The Topology must be a derived spatially decomposed type (not Rep_Topology).

 * \sa rtt_mesh::Comm_Patterns and rtt_mesh::Topology (rtt_mesh::DD_Topology)
 *
 * \param[in] pattern active Comm_Patterns object, the patterns must be set
 * \param[in] local_data field of local data on processor that will be sent to other processors
 *               according to the Comm_Patterns. Must implement operator[] with the cell-1
 *               convention.
 * \param[in,out] ghost_data empty ghost cell field that will be filled. Must have size equal to the
 *               number of ghost nodes on the local processor. Must implement operator[] with the
 *               cell-1 convention.
 * \param[in] vector_size  number of DoFs per cell.
 *
 *  NOTE: this function assumes 1D array, size of cell*vector_size.
 *        Default is scalar (vector_size=1)
 */
template <typename FT, typename B_IT>
void Parallel_Data_Operator::gather_ghost_cell_data(SP_Comm_Patterns const &pattern,
                                                    const FT &local_data, B_IT ghost_data,
                                                    const size_t vector_size) const {
  Require(pattern);
  Require(local_data.size() > 0);
  // This is a different tag than then parallel operation of serial values though we don't think
  // these function calls will ever conflict
  static const int comm_tag = 606;

  // Extract the type of the data being received. Data being sent will be converted to this type.
  using Vector_T = typename std::iterator_traits<B_IT>::value_type;
  using Vector_Iter_T = typename Vector_T::iterator;
  using Recv_T = typename std::iterator_traits<Vector_Iter_T>::value_type;
  using Data = Data_Request<Recv_T>;

  Check(*pattern);

  const size_t num_incoming = pattern->get_num_recv_procs();

  // Post receives to the data we need.
  // ----------------------------------

  // Make a vector of data that will be received from each processor: We will stack cell wise vector
  // data to a 1d vector to send and receive the data and we will then repopulate it into the proper
  // vector structure
  std::vector<Data> incoming_data(num_incoming, Data());

  auto const recv_begin = pattern->get_recv_begin();
  auto const recv_end = pattern->get_recv_end();

  auto data_request = incoming_data.begin();

  for (auto itor = recv_begin; itor != recv_end; ++itor, ++data_request) {
    // Determine the processor to receive from and the size of the data.
    int const processor = itor->first;
    size_t const data_size = itor->second.size();

    data_request->allocate(data_size * vector_size);
    data_request->post(processor, comm_tag);
  }
  Check(data_request == incoming_data.end());

  // Send data we have that other processors need.
  // ---------------------------------------------

  auto const send_begin = pattern->get_send_begin();
  auto const send_end = pattern->get_send_end();

  for (auto itor = send_begin; itor != send_end; itor++) {
    // Determine the processor we are sending to and the size of the data.
    int const processor = itor->first;
    size_t const data_size = itor->second.size();

    Data_Send<Recv_T> sender(data_size * vector_size);

    // Pack the data to send into the sender's array.
    Recv_T *data = sender.get_data();
    size_t k = 0;
    for (size_t i = 0; i < data_size; i++) {
      for (size_t j = 0; j < vector_size; j++) {
        // Get the local cell index and assign the data
        int const local_cell = itor->second[i];
        Check(local_cell > 0);
        Check(local_data[local_cell - 1].size() == vector_size);
        data[k] = local_data[local_cell - 1][j];
        k++;
      }
    }
    Check(k == data_size * vector_size);
    sender.send(processor, comm_tag);
  }

  // Receive the data and write out ghost_data field
  // --------------------------------------------

  size_t finished = 0;
  while (finished < num_incoming) {
    // Loop through the C4 requests for each processor that we communicate with and see if the data
    // is here

    data_request = incoming_data.begin();
    for (auto itor = recv_begin; itor != recv_end; ++itor, ++data_request) {
      size_t const data_size = itor->second.size();

      // If the message is complete...
      if (data_request->complete()) {
        Require(!data_request->inuse());

        finished++;

        // Write the data to the bc_data field
        Recv_T *data = data_request->get_data();
        size_t k = 0;
        for (size_t c = 0; c < data_size; ++c) {
          // Determine the boundary cell and write the data
          int const ghost_cell = itor->second[c];
          for (size_t j = 0; j < vector_size; j++) {
            ghost_data[ghost_cell - 1][j] = data[k];
            k++;
          }
        }
        Check(k == data_size * vector_size);

        data_request->release();
      }
    }
    Check(data_request == incoming_data.end());
  }

  Ensure(finished == num_incoming);

  // sync everything
  rtt_c4::global_barrier();
}

//------------------------------------------------------------------------------------------------//
// GHOST CELL FLAT VECTOR PARALLEL OPERATIONS (GATHER/SCATTER TYPE)
//------------------------------------------------------------------------------------------------//
/*!
 * \brief Fill boundary cell vector data fields in a spatially decomposed topology.
 *
 * This function utilizes the rtt_mesh::Comm_Patterns class to "gather" data into ghost cell
 * fields.  What this entails is the following:
 *
 * \arg using the Comm_Patterns, each processor posts requests to processors it needs data from;

 * \arg using the Comm_Patterns and local data fields, each processor sends data out to processors
 * that need data;

 * \arg each processor receives the data for its boundary cells and writes it into a local boundary
 * cell field

 * The Topology must be a derived spatially decomposed type (not Rep_Topology).

 * \sa rtt_mesh::Comm_Patterns and rtt_mesh::Topology (rtt_mesh::DD_Topology)
 *
 * \param[in] pattern active Comm_Patterns object, the patterns must be set
 * \param[in] local_data field of local data on processor that will be sent to other processors
 *               according to the Comm_Patterns. Must implement operator[] with the cell-1
 *               convention.
 * \param[in,out] ghost_data empty ghost cell field that will be filled. Must have size equal to the
 *               number of ghost nodes on the local processor. Must implement operator[] with the
 *               cell-1 convention.
 * \param[in] vector_size  number of DoFs per cell.
 *
 *  NOTE: this function assumes 1D array, size of cell*vector_size.
 *        Default is scalar (vector_size=1)
 */
template <typename FT, typename B_IT>
void Parallel_Data_Operator::gather_fv_ghost_cell_data(SP_Comm_Patterns const &pattern,
                                                       const FT &local_data, B_IT ghost_data,
                                                       const size_t vector_size) const {
  Require(pattern);
  Require(local_data.size() > 0);
  // This is a different tag than then parallel operation of serial values though we don't think
  // these function calls will ever conflict
  static const int comm_tag = 606;

  // Extract the type of the data being received. Data being sent will be converted to this type.
  using Recv_T = typename std::iterator_traits<B_IT>::value_type;
  using Data = Data_Request<Recv_T>;

  Check(*pattern);

  const size_t num_incoming = pattern->get_num_recv_procs();

  // Post receives to the data we need.
  // ----------------------------------

  // Make a vector of data that will be received from each processor: We will stack cell wise vector
  // data to a 1d vector to send and receive the data and we will then repopulate it into the proper
  // vector structure
  std::vector<Data> incoming_data(num_incoming, Data());

  auto const recv_begin = pattern->get_recv_begin();
  auto const recv_end = pattern->get_recv_end();

  auto data_request = incoming_data.begin();

  for (auto itor = recv_begin; itor != recv_end; ++itor, ++data_request) {
    // Determine the processor to receive from and the size of the data.
    int const processor = itor->first;
    size_t const data_size = itor->second.size();

    data_request->allocate(data_size * vector_size);
    data_request->post(processor, comm_tag);
  }
  Check(data_request == incoming_data.end());

  // Send data we have that other processors need.
  // ---------------------------------------------

  auto const send_begin = pattern->get_send_begin();
  auto const send_end = pattern->get_send_end();

  for (auto itor = send_begin; itor != send_end; itor++) {
    // Determine the processor we are sending to and the size of the data.
    int const processor = itor->first;
    size_t const data_size = itor->second.size();

    Data_Send<Recv_T> sender(data_size * vector_size);

    // Pack the data to send into the sender's array.
    Recv_T *data = sender.get_data();
    size_t k = 0;
    for (size_t i = 0; i < data_size; i++) {
      for (size_t j = 0; j < vector_size; j++) {
        // Get the local cell index and assign the data
        int const local_cell = itor->second[i];
        Check(local_cell > 0);
        data[k] = local_data[(local_cell - 1) * vector_size + j];
        k++;
      }
    }
    Check(k == data_size * vector_size);
    sender.send(processor, comm_tag);
  }

  // Receive the data and write out ghost_data field
  // --------------------------------------------

  size_t finished = 0;
  while (finished < num_incoming) {
    // Loop through the C4 requests for each processor that we communicate with and see if the data
    // is here

    data_request = incoming_data.begin();
    for (auto itor = recv_begin; itor != recv_end; ++itor, ++data_request) {
      size_t const data_size = itor->second.size();

      // If the message is complete...
      if (data_request->complete()) {
        Require(!data_request->inuse());

        finished++;

        // Write the data to the bc_data field
        Recv_T *data = data_request->get_data();
        size_t k = 0;
        for (size_t c = 0; c < data_size; ++c) {
          // Determine the boundary cell and write the data
          int const ghost_cell = itor->second[c];
          for (size_t j = 0; j < vector_size; j++) {
            ghost_data[(ghost_cell - 1) * vector_size + j] = data[k];
            k++;
          }
        }
        Check(k == data_size * vector_size);

        data_request->release();
      }
    }
    Check(data_request == incoming_data.end());
  }

  Ensure(finished == num_incoming);

  // sync everything
  rtt_c4::global_barrier();
}

//------------------------------------------------------------------------------------------------//
/*!

 * \brief Fill boundary cell data fields in a spatially decomposed topology.

 * This function utilizes the rtt_mc::Comm_Patterns class to "gather" data into boundary cell
 * fields.  What this entails is the following:

 * \arg using the Comm_Patterns, each processor posts requests to processors it needs data from;

 * \arg using the Comm_Patterns and local data fields, each processor sends data out to processors
 *      that need data;

 * \arg each processor receives the data for its boundary cells and writes it into a local boundary
 *      cell field

 * The Topology must be a derived spatially decomposed type (not Rep_Topology).

 * \sa rtt_mc::Comm_Patterns and rtt_mc::Topology (rtt_mc::DD_Topology)
 *
 * \param pattern active Comm_Patterns object, the patterns must be set
 *
 * \param local_data field of local data on processor that will be sent to other processors
 *               according to the Comm_Patterns. Must implement operator[] with the cell-1
 *               convention.
 *
 * \param bc_data empty boundary cell field that will be filled. Must have size equal to the number
 *               of boundary nodes on the local processor. Must implement operator[] with the cell-1
 *               convention.
 */
template <typename FT, typename B_IT>
void Parallel_Data_Operator::gather_ghost_cell_data(SP_Comm_Patterns const &pattern,
                                                    const FT &local_data, B_IT bc_data) const {
  static const int comm_tag = 605;

  // Extract the type of the data being received. Data being sent will be converted to this type.
  using Recv_T = typename std::iterator_traits<B_IT>::value_type;
  using Data = Data_Request<Recv_T>;

  Check(*pattern);

  const size_t num_incoming = pattern->get_num_recv_procs();

  // Post receives to the data we need.
  // ----------------------------------

  // Make a vector of data that will be received from each processor:
  std::vector<Data> incoming_data(num_incoming, Data());

  auto const recv_begin = pattern->get_recv_begin();
  auto const recv_end = pattern->get_recv_end();

  auto data_request = incoming_data.begin();

  for (auto itor = recv_begin; itor != recv_end; ++itor, ++data_request) {
    // Determine the processor to receive from and the size of the data.
    int const processor = itor->first;
    size_t const data_size = itor->second.size();

    data_request->allocate(data_size);
    data_request->post(processor, comm_tag);
  }
  Check(data_request == incoming_data.end());

  // Send data we have that other processors need.
  // ---------------------------------------------

  auto const send_begin = pattern->get_send_begin();
  auto const send_end = pattern->get_send_end();

  for (auto itor = send_begin; itor != send_end; itor++) {
    // Determine the processor we are sending to and the size of the data.
    int const processor = itor->first;
    size_t const data_size = itor->second.size();

    Data_Send<Recv_T> sender(data_size);

    // Pack the data to send into the sender's array.
    Recv_T *data = sender.get_data();
    for (size_t i = 0; i < data_size; i++) {
      // Get the local cell index and assign the data
      int const local_cell = itor->second[i];
      Check(local_cell > 0);
      data[i] = local_data[local_cell - 1];
    }

    sender.send(processor, comm_tag);
  }

  // Receive the data and write out bc_data field
  // --------------------------------------------

  size_t finished = 0;
  while (finished < num_incoming) {
    // Loop through the C4 requests for each processor that we communicate with and see if the data
    // is here

    data_request = incoming_data.begin();
    for (auto itor = recv_begin; itor != recv_end; ++itor, ++data_request) {
      size_t const data_size = itor->second.size();

      // If the message is complete...
      if (data_request->complete()) {
        Require(!data_request->inuse());

        finished++;

        // Write the data to the bc_data field
        Recv_T *data = data_request->get_data();
        for (size_t c = 0; c < data_size; ++c) {
          // Determine the boundary cell and write the data
          int const boundary_cell = itor->second[c];
          bc_data[boundary_cell - 1] = data[c];
        }

        data_request->release();
      }
    }
    Check(data_request == incoming_data.end());
  }

  Ensure(finished == num_incoming);

  // sync everything
  rtt_c4::global_barrier();
}

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Fill boundary cell data fields in a spatially decomposed topology.
 *
 * This function utilizes the rtt_mesh::Comm_Patterns class to "gather" data into boundary cell
 * fields.  What this entails is the following:
 *
 * \arg using the Comm_Patterns, each processor posts requests to processors it needs data from;
 *
 * \arg using the Comm_Patterns and local data fields, each processor sends data out to processors
 *      that need data;
 *
 * \arg each processor receives the data for its boundary cells and writes it into a local boundary
 *      cell field
 *
 * The Topology must be a derived spatially decomposed type (not Rep_Topology).
 *
 * \sa rtt_mesh::Comm_Patterns and rtt_mesh::Topology (rtt_mesh::DD_Topology)
 *
 * \param pattern active Comm_Patterns object, the patterns must be set
 * \param local_data vector of local data on processor that will be sent to other processors
 *               according to the Comm_Patterns
 * \param bc_data empty boundary cell field that will be filled
 */
template <typename T>
void Parallel_Data_Operator::gather_ghost_cell_data(SP_Comm_Patterns pattern,
                                                    const std::vector<T> &local_data,
                                                    std::vector<T> &bc_data) {
  Require(bc_data.empty());
  Require(topology->num_local_cells() < std::numeric_limits<int>::max());
  Require(local_data.size() == static_cast<size_t>(topology->num_local_cells()));
  Require(*pattern);

  bc_data.resize(topology->num_ghost_cells());

  gather_ghost_cell_data(pattern, local_data, bc_data.begin());
}

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Fill boundary cell vector data fields in a spatially decomposed topology.
 *
 * This function utilizes the rtt_mesh::Comm_Patterns class to "gather" data into boundary cell
 * fields.  What this entails is the following:
 *
 * \arg using the Comm_Patterns, each processor posts requests to processors it needs data from;
 *
 * \arg using the Comm_Patterns and local data fields, each processor sends data out to processors
 *      that need data;
 *
 * \arg each processor receives the data for its boundary cells and writes it into a local boundary
 *      cell field
 *
 * The Topology must be a derived spatially decomposed type (not Rep_Topology).
 *
 * \sa rtt_mesh::Comm_Patterns and rtt_mesh::Topology (rtt_mesh::DD_Topology)
 *
 * \param[in] pattern active Comm_Patterns object, the patterns must be set
 * \param[in] local_data vector of local data on processor that will be sent to other processors
 *               according to the Comm_Patterns
 * \param[in,out] bc_data empty boundary cell field that will be filled
 * \param[in] vector_size
 * \bug Need to document parameter vector_size
 */
template <typename T>
void Parallel_Data_Operator::gather_ghost_cell_data(SP_Comm_Patterns pattern,
                                                    const std::vector<std::vector<T>> &local_data,
                                                    std::vector<std::vector<T>> &bc_data,
                                                    size_t const vector_size) {
  Require(bc_data.empty());
  Require(topology->num_local_cells() < std::numeric_limits<int>::max());
  Require(local_data.size() == static_cast<size_t>(topology->num_local_cells()));
  Require(*pattern);

  bc_data.resize(topology->num_ghost_cells(), std::vector<T>(vector_size));

  gather_ghost_cell_data(pattern, local_data, bc_data.begin(), vector_size);
}

//------------------------------------------------------------------------------------------------//
// GHOST CELL RAGGED PARALLEL OPERATIONS (GATHER/SCATTER TYPE)
//------------------------------------------------------------------------------------------------//
/*!
 * \brief Fill boundary cell vector data fields in a spatially decomposed topology.
 *
 * This function utilizes the rtt_mesh::Comm_Patterns class to "gather" data into ghost cell
 * fields.  What this entails is the following:
 *
 * \arg using the Comm_Patterns, each processor posts requests to processors it needs data from;

 * \arg using the Comm_Patterns and local data fields, each processor sends data out to processors
 * that need data;

 * \arg each processor receives the data for its boundary cells and writes it into a local boundary
 * cell field

 * The Topology must be a derived spatially decomposed type (not Rep_Topology).

 * \sa rtt_mesh::Comm_Patterns and rtt_mesh::Topology (rtt_mesh::DD_Topology)
 *
 * \param[in] pattern active Comm_Patterns object, the patterns must be set
 * \param[in] local_data field of local data on processor that will be sent to other processors
 *               according to the Comm_Patterns. Must implement operator[] with the cell-1
 *               convention.
 * \param[in,out] ghost_data empty ghost cell field that will be filled. Must have size equal to the
 *               number of ghost nodes on the local processor for each receiving rank. Must
 *               implement operator[] with the cell-1 convention.
 * \param[in] local_vector_size vector of number of DoFs per local cell
 * \param[in] ghost_vector_size nested vector of number of DoFs per ghost cell per off-processor
              rank.
 *
 *  NOTE: this function assumes 1D array, size of cell*vector_size.
 *        Default is scalar (vector_size=1)
 */
template <typename T>
void Parallel_Data_Operator::gather_ghost_cell_data(
    SP_Comm_Patterns const &pattern, const std::vector<std::vector<T>> &local_data,
    std::vector<std::vector<T>> &ghost_data, const std::vector<size_t> &local_vector_size,
    const std::vector<size_t> &ghost_vector_size) const {

  Require(pattern);
  Require(ghost_data.empty());
  Require(local_data.size() > 0);
  // This is a different tag than the parallel operation of serial/fixed vector values though we
  // don't think these function calls will ever conflict
  static const int comm_tag = 607;

  // Allocate ghost data memory
  ghost_data.resize(topology->num_ghost_cells());
  for (size_t cell = 0; cell < topology->num_ghost_cells(); cell++) {
    ghost_data[cell] = std::vector<T>(ghost_vector_size[cell]);
  }

  // Extract the type of the data being received. Data being sent will be converted to this type.
  using Recv_T = typename std::decay_t<decltype(*(*ghost_data.begin()).begin())>;
  using Data = Data_Request<Recv_T>;

  Check(*pattern);

  const size_t num_incoming = pattern->get_num_recv_procs();

  // Post receives to the data we need.
  // ----------------------------------

  // Make a vector of data that will be received from each processor: We will stack cell wise vector
  // data to a 1d vector to send and receive the data and we will then repopulate it into the proper
  // vector structure
  std::vector<Data> incoming_data(num_incoming, Data());

  auto const recv_begin = pattern->get_recv_begin();
  auto const recv_end = pattern->get_recv_end();

  auto data_request = incoming_data.begin();

  for (auto itor = recv_begin; itor != recv_end; ++itor, ++data_request) {
    // Determine the processor to receive from and the size of the data.
    int const processor = itor->first;

    // Processor ordering of comm_patterns map must match provided ghost_vector_size.
    size_t data_size = 0;
    for (auto ghost_cell : itor->second) {
      data_size += ghost_vector_size[ghost_cell - 1];
    }

    data_request->allocate(data_size);
    data_request->post(processor, comm_tag);
  }
  Check(data_request == incoming_data.end());

  // Send data we have that other processors need.
  // ---------------------------------------------

  auto const send_begin = pattern->get_send_begin();
  auto const send_end = pattern->get_send_end();

  for (auto itor = send_begin; itor != send_end; itor++) {
    // Determine the processor we are sending to and the size of the data.
    int const processor = itor->first;
    size_t data_size = 0;
    for (auto cell : itor->second) {
      data_size += local_vector_size[cell - 1];
    }
    Require(data_size > 0);

    Data_Send<Recv_T> sender(data_size);

    // Pack the data to send into the sender's array.
    Recv_T *data = sender.get_data();
    size_t k = 0;
    for (size_t cell = 0; cell < itor->second.size(); cell++) {
      for (size_t j = 0; j < local_vector_size[cell]; j++) {
        int const local_cell = itor->second[cell];
        Check(local_cell > 0);
        Check(local_data[local_cell - 1].size() == local_vector_size[cell]);
        data[k] = local_data[local_cell - 1][j];
        k++;
      }
    }
    sender.send(processor, comm_tag);
  }

  // Receive the data and write out ghost_data field
  // --------------------------------------------

  size_t finished = 0;
  while (finished < num_incoming) {
    // Loop through the C4 requests for each processor that we communicate with and see if the data
    // is here

    data_request = incoming_data.begin();
    for (auto itor = recv_begin; itor != recv_end; ++itor, ++data_request) {

      // If the message is complete...
      if (data_request->complete()) {
        Require(!data_request->inuse());

        finished++;

        // Write the data to the bc_data field
        Recv_T *data = data_request->get_data();
        size_t k = 0;
        for (auto ghost_cell : itor->second) {
          // Determine the boundary cell and write the data
          for (size_t j = 0; j < ghost_vector_size[ghost_cell - 1]; j++) {
            ghost_data[ghost_cell - 1][j] = data[k];
            k++;
          }
        }

        data_request->release();
      }
    }
    Check(data_request == incoming_data.end());
  }

  Ensure(finished == num_incoming);

  // sync everything
  rtt_c4::global_barrier();
}

//------------------------------------------------------------------------------------------------//
// GHOST CELL RAGGED VECTOR PARALLEL OPERATIONS (GATHER/SCATTER TYPE)
//------------------------------------------------------------------------------------------------//
/*!
 * \brief Fill boundary cell vector data fields in a spatially decomposed topology.
 *
 * This function utilizes the rtt_mesh::Comm_Patterns class to "gather" data into ghost cell
 * fields.  What this entails is the following:
 *
 * \arg using the Comm_Patterns, each processor posts requests to processors it needs data from;

 * \arg using the Comm_Patterns and local data fields, each processor sends data out to processors
 * that need data;

 * \arg each processor receives the data for its boundary cells and writes it into a local boundary
 * cell field

 * The Topology must be a derived spatially decomposed type (not Rep_Topology).

 * \sa rtt_mesh::Comm_Patterns and rtt_mesh::Topology (rtt_mesh::DD_Topology)
 *
 * \param[in] pattern active Comm_Patterns object, the patterns must be set
 * \param[in] local_data field of local data on processor that will be sent to other processors
 *               according to the Comm_Patterns. Must implement operator[] with the cell-1
 *               convention.
 * \param[in,out] ghost_data empty ghost cell field that will be filled. Must have size equal to the
 *               number of ghost nodes on the local processor for each receiving rank. Must
 *               implement operator[] with the cell-1 convention.
 * \param[in] local_vector_size vector of number of DoFs per local cell
 * \param[in] ghost_vector_size nested vector of number of DoFs per ghost cell per off-processor
              rank.
 * \param[in] Nv constant length vector on top of per-cell multiplicity
 *
 *  NOTE: this function assumes 1D array, size of cell*vector_size.
 *        Default is scalar (vector_size=1)
 */
template <typename T>
void Parallel_Data_Operator::gather_ghost_cell_data(
    SP_Comm_Patterns const &pattern, const std::vector<std::vector<std::vector<T>>> &local_data,
    std::vector<std::vector<std::vector<T>>> &ghost_data,
    const std::vector<size_t> &local_vector_size, const std::vector<size_t> &ghost_vector_size,
    const int Nv) const {
  Require(pattern);
  Require(ghost_data.empty());
  Require(local_data.size() > 0);
  // This is a different tag than the parallel operation of serial/fixed vector values though we
  // don't think these function calls will ever conflict
  static const int comm_tag = 608;

  // Allocate ghost data memory
  ghost_data.resize(topology->num_ghost_cells());
  for (size_t cell = 0; cell < topology->num_ghost_cells(); cell++) {
    ghost_data[cell] = std::vector<std::vector<T>>(ghost_vector_size[cell]);
    for (size_t m = 0; m < ghost_vector_size[cell]; m++) {
      ghost_data[cell][m] = std::vector<T>(Nv);
    }
  }

  // Extract the type of the data being received. Data being sent will be converted to this type.
  using Recv_T = typename std::decay_t<decltype(*(*(*ghost_data.begin()).begin()).begin())>;
  using Data = Data_Request<Recv_T>;

  Check(*pattern);

  const size_t num_incoming = pattern->get_num_recv_procs();

  // Post receives to the data we need.
  // ----------------------------------

  // Make a vector of data that will be received from each processor: We will stack cell wise vector
  // data to a 1d vector to send and receive the data and we will then repopulate it into the proper
  // vector structure
  std::vector<Data> incoming_data(num_incoming, Data());

  auto const recv_begin = pattern->get_recv_begin();
  auto const recv_end = pattern->get_recv_end();

  auto data_request = incoming_data.begin();

  for (auto itor = recv_begin; itor != recv_end; ++itor, ++data_request) {
    // Determine the processor to receive from and the size of the data.
    int const processor = itor->first;

    // Processor ordering of comm_patterns map must match provided ghost_vector_size.
    //Require(itor->second.size() == ghost_vector_size.size());
    size_t data_size = 0;
    for (auto ghost_cell : itor->second) {
      data_size += ghost_vector_size[ghost_cell - 1];
    }
    data_size *= Nv;

    data_request->allocate(data_size);
    data_request->post(processor, comm_tag);
  }
  Check(data_request == incoming_data.end());

  // Send data we have that other processors need.
  // ---------------------------------------------

  auto const send_begin = pattern->get_send_begin();
  auto const send_end = pattern->get_send_end();

  for (auto itor = send_begin; itor != send_end; itor++) {
    // Determine the processor we are sending to and the size of the data.
    int const processor = itor->first;
    size_t data_size = 0;
    for (auto cell : itor->second) {
      data_size += local_vector_size[cell - 1];
    }
    data_size *= Nv;
    Require(data_size > 0);

    Data_Send<Recv_T> sender(data_size);

    // Pack the data to send into the sender's array.
    Recv_T *data = sender.get_data();
    size_t k = 0;
    for (size_t cell = 0; cell < itor->second.size(); cell++) {
      for (size_t j = 0; j < local_vector_size[cell]; j++) {
        for (auto n = 0; n < Nv; n++) {
          int const local_cell = itor->second[cell];
          Check(local_cell > 0);
          Check(local_data[local_cell - 1].size() == local_vector_size[cell]);
          data[k] = local_data[local_cell - 1][j][n];
          k++;
        }
      }
    }
    sender.send(processor, comm_tag);
  }

  // Receive the data and write out ghost_data field
  // --------------------------------------------

  size_t finished = 0;
  while (finished < num_incoming) {
    // Loop through the C4 requests for each processor that we communicate with and see if the data
    // is here

    data_request = incoming_data.begin();
    for (auto itor = recv_begin; itor != recv_end; ++itor, ++data_request) {

      // If the message is complete...
      if (data_request->complete()) {
        Require(!data_request->inuse());

        finished++;

        // Write the data to the bc_data field
        Recv_T *data = data_request->get_data();
        size_t k = 0;
        for (auto ghost_cell : itor->second) {
          // Determine the boundary cell and write the data
          for (size_t j = 0; j < ghost_vector_size[ghost_cell - 1]; j++) {
            for (int n = 0; n < Nv; n++) {
              ghost_data[ghost_cell - 1][j][n] = data[k];
              k++;
            }
          }
        }

        data_request->release();
      }
    }
    Check(data_request == incoming_data.end());
  }

  Ensure(finished == num_incoming);

  // sync everything
  rtt_c4::global_barrier();
}

////////////////////////////////////////////////////////////
// ghost cell data exchange
// scalar case
// local data will be resized and append the ghost values..
template <typename T>
void Parallel_Data_Operator::fill_ghost_cell_data(SP_Comm_Patterns pattern,
                                                  std::vector<T> &local_data,
                                                  const size_t vector_size) {
  size_t const nowned = topology->num_local_cells();
  size_t const numGhosts = topology->num_ghost_cells();
  if (numGhosts == 0)
    return;

  std::vector<T> ghost_data(numGhosts * vector_size);

  gather_fv_ghost_cell_data(pattern, local_data, ghost_data.begin(), vector_size);
  size_t index(0);
  // resize local_data..
  local_data.resize((nowned + numGhosts) * vector_size);

  for (size_t i = nowned; i < nowned + numGhosts; ++i) {
    for (size_t j = 0; j < vector_size; ++j)
      local_data[i * vector_size + j] = ghost_data[index++];
  }
}
} // end namespace rtt_mesh

#endif // mesh_Parallel_Data_Operator_i_hh

//------------------------------------------------------------------------------------------------//
// end of mesh/Parallel_Data_Operator.i.hh
//------------------------------------------------------------------------------------------------//
