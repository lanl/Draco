//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   mesh/CAMR_Mesh.t.hh
 * \author H. Park
 * \date   May 16 2023
 * \brief  Member definitions of class CAMR_Mesh
 * \note   Copyright (C) 2023-2024 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#ifndef rtt_mesh_CAMR_Mesh_t_hh
#define rtt_mesh_CAMR_Mesh_t_hh

#include "CAMR_Mesh.hh"
#include "viz/Ensight_Translator.hh"
#include <bitset>
#include <limits>

namespace rtt_mesh {

//------------------------------------------------------------------------------------------------//
/**
 * \brief Return the width of the given level of the mesh
 * \arg level The level of the mesh. Defaults to zero.
 */
template <int D> double CAMR_Mesh<D>::get_width(size_t level) const {
  //Check(level >= 0);
  Check(level < number_of_levels);
  return widths[level];
}

//------------------------------------------------------------------------------------------------//
/**
 * \brief Return the width of the given cell in the mesh
 *
 * \arg cell The cell in the mesh.
 * \return The uniform width of the cell.
 */
template <int D> double CAMR_Mesh<D>::get_cell_width(size_t cell) const {
  Require(cell >= 1 && cell <= levels.size());
  return get_width(get_level(cell));
}

//------------------------------------------------------------------------------------------------//
//! Set values for the data members related to the widths of the various levels of the mesh.
template <int D> void CAMR_Mesh<D>::build_widths() {
  maximum_level = *std::max_element(levels.begin(), levels.end());
  number_of_levels = maximum_level + 1;

  // Require (maximum_level >= 0);
  // Compute the contents of the widths array.
  widths.resize(number_of_levels);

  double power_of_two = 1.0;
  for (size_t level = 0; level <= maximum_level; ++level) {
    widths[level] = width * power_of_two;
    power_of_two /= 2.0;
  }
}

// Suppress the warning: 'conditional expression is constant'. This can be removed once we adopt C++
// and the syntax is changed to 'if constexpr'
#if defined(MSVC)
#pragma warning(push)
#pragma warning(disable : 4127)
#endif

//------------------------------------------------------------------------------------------------//
//! Compute the total volume covered by this mesh.
template <int D> void CAMR_Mesh<D>::compute_total_volume() {
  Check(active_size < INT32_MAX);
  for (size_t cell = 1; cell <= active_size; ++cell) {
    double const cell_width = get_cell_width(static_cast<int>(cell));
    double cell_volume = cell_width * cell_width;
    //! \bug This should be 'if constexpr' (C++17)
    if (D == 3)
      cell_volume *= cell_width;
    local_volume += cell_volume;
  }
}

#if defined(MSVC)
#pragma warning(pop)
#endif

//------------------------------------------------------------------------------------------------//
/**
 * \brief Return the level of a cell
 * \arg The cell
 * \return The level
 */
template <int D> size_t CAMR_Mesh<D>::get_level(size_t cell) const {
  Require(cell >= 1 && cell <= levels.size());
  Ensure(static_cast<uint32_t>(levels[cell - 1]) < UINT32_MAX);
  return static_cast<size_t>(levels[cell - 1]);
}

//------------------------------------------------------------------------------------------------//
/**
 * \brief Return a vector containing the canonical coordinates for a given cell
 *
 * \arg cell  The cell in question
 *
 * \return A std::vector<double> size D, containing the x,y[,z] coordinates of the vertex of the
 *         cell with the lowest value in all widths. (This is the canonical vertex of the cell).
 */
template <int D> typename CAMR_Mesh<D>::sf_double CAMR_Mesh<D>::get_vertex(size_t cell) const {
  Require(cell >= 1 && cell <= layout->get_num_total_cells());
  auto data = vertex_data(cell);

  sf_double vertex(data, data + D);

  Ensure(vertex.size() == D);
  return vertex;
}
//------------------------------------------------------------------------------------------------//
/*!
 * \brief Return the minimum cell boundary along dimension d in a cell.
 *
 * \param d mesh dimension (1,2,3) for (x,y,z)
 * \param cell cell index [1,N]
 */
template <int D> double CAMR_Mesh<D>::min(size_t d, size_t cell) const {
  Require(d > 0 && d <= D);
  Require(cell > 0 && static_cast<uint32_t>(cell) <= layout->get_num_total_cells());

  Check(cell < INT32_MAX);
  sf_double const cell_vertex = get_vertex(static_cast<int>(cell));

  double minimum(cell_vertex[d - 1]);
  // return minimum dimension
  return minimum;
}

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Return the maximum cell boundary along dimension d in a cell.
 *
 * \param d mesh dimension (1,2,3) for (x,y,z)
 * \param cell cell index [1,N]
 */
template <int D> double CAMR_Mesh<D>::max(size_t d, size_t cell) const {
  Require(d > 0 && d <= D);
  Require(cell > 0 && static_cast<uint32_t>(cell) <= layout->get_num_total_cells());
  Check(cell < INT32_MAX);
  sf_double const cell_vertex = get_vertex(static_cast<int>(cell));
  double const delta = get_cell_width(static_cast<int>(cell));

  double maximum(cell_vertex[d - 1] + delta);
  // return maximum dimension
  return maximum;
}

//------------------------------------------------------------------------------------------------//
/**
 * \brief Return true if the given point is within the given cell based on a soft check.
 *
 * \arg cell The cell
 * \arg position The point
 * \return boolean, True if the point is within the cell.
 */
template <int D> bool CAMR_Mesh<D>::in_cell(size_t cell, const U3 &position) const {
  using rtt_dsxx::soft_equiv;
  Require(cell > 0 && static_cast<uint32_t>(cell) <= layout->get_num_owned_cells());
  sf_double const cell_vertex = get_vertex(cell);
  double const delta = get_cell_width(cell);

  // high and low sides
  double high = 0.0;
  double low = 0.0;

  for (int d = 0; d < D; ++d) {
    // assign low and high sides for this dimension
    high = cell_vertex[d] + delta;
    low = cell_vertex[d];

    // low side check
    if (position[d] < low && !soft_equiv(position[d], low))
      return false;

    // high side
    if (position[d] > high && !soft_equiv(position[d], high))
      return false;
  }

  return true;
}

//------------------------------------------------------------------------------------------------//
/**
 * \brief Return true if the given point is strictly within the given cell based on a hard check.
 *
 * \arg cell The cell
 * \arg position The point
 * \return boolean, True if the point is within the cell.
 */
template <int D> bool CAMR_Mesh<D>::hard_in_cell(size_t const cell, const U3 &position) const {
  Require(cell > 0 && static_cast<uint32_t>(cell) <= layout->get_num_owned_cells());
  sf_double const cell_vertex = get_vertex(cell);
  double const delta = get_cell_width(cell);

  for (uint32_t dimension = 1; dimension <= D; ++dimension) {
    uint32_t dimension_index = dimension - 1;
    if (position[dimension_index] < cell_vertex[dimension_index] ||
        position[dimension_index] > cell_vertex[dimension_index] + delta)
      return false;
  }
  return true;
}

//------------------------------------------------------------------------------------------------//
/**
 * \brief Return true if the given point is on the "inside" side of the given cell with respect to
 *        the given face.
 *
 * \arg cell The cell
 * \arg face The face
 * \arg position The point
 * \return boolean, true if point is on the "inside" side of the face side of cell.
 */
template <int D>
bool CAMR_Mesh<D>::in_cell(size_t const cell, uint32_t face, const U3 &position) const {
  Require(cell > 0 && static_cast<uint32_t>(cell) <= layout->get_num_owned_cells());
  Require(face >= 1 && face <= layout->num_coarse_faces());

  //  Direction<D> d(face);
  // check if face is positive, (ordered in -x,+x, -y,+y, -z,+z)
  bool const positive = static_cast<bool>((face - 1) % 2); //d.get_positive();
  // FIXME: I don't think this is correct??
  int const dimension_index = (face - 1) % 2; // d.get_axis();

  double const delta = get_cell_width(cell);
  double const pos = position[dimension_index];
  double const vertex = get_vertex(cell)[dimension_index];

  return (positive && (pos < vertex + delta)) || (!positive && (pos > vertex));
}

//------------------------------------------------------------------------------------------------//
// Interface Functions
//------------------------------------------------------------------------------------------------//
/**
 * \brief Return the cell across a face from a given cell. Use the position argument to resolve the
 *        correct cell if more than one is present.
 *
 * \arg cell  The cell in question
 * \arg face  The face in question
 * \arg position  A position vector for a point on (or very near) the face.
 */
template <int D>
int CAMR_Mesh<D>::next_cell(const int cell, const uint32_t face, const U3 &position) const {
  Require(cell > 0 && static_cast<uint32_t>(cell) <= layout->get_num_total_cells());
  Require(face >= 1 && face <= layout->num_coarse_faces());

  if (layout->num_cells_across(cell, face) == 1) {
    // all AMR faces are the same, just return the first one
    Check((*layout)(cell, face, 1) < INT32_MAX);
    return (*layout)(cell, face, 1);
  }

  /* The neighboring face is refined. To determine the next cell, we extract the D-1 coordinates of
   * the point that are on the face being crossed and determine if they are above or below the
   * midpoint of the cell. From this information, we select the next cell among the neighbors.
   *
   * This is where we implicitly use the assumption that the refined neighbors of a cell are listed
   * in x,y,z order.
   *
   * XXX Try to move this assumption into layout itself. Have layout return the next cell based on
   * the "above/below midpoint" coordinates. (Called "neighbor_coords" in the code below).
   */

  std::pair<int, bool> const face_data = layout->face_axis_sign(face);
  int const &dimension = face_data.first;

  sf_double const cell_coord(get_vertex(cell));
  double delta = get_cell_width(cell);

  std::array<int, 2> neighbor_coords = {0, 0};
  auto insert_coord = neighbor_coords.begin();

  for (int other_dimension = 1; other_dimension <= D; ++other_dimension) {
    if (other_dimension == dimension)
      continue;

    const double relative_position =
        (position[other_dimension - 1] - cell_coord[other_dimension - 1]) / delta;

    Require(relative_position >= 0.0);
    Require(relative_position <= 1.0);

    *insert_coord++ = static_cast<int>(relative_position > 0.5);
  }

  int const neighbor_index = neighbor_coords[0] + 2 * neighbor_coords[1];

  Check((*layout)(cell, face)[neighbor_index] < INT32_MAX);
  auto next_cell = (*layout)(cell, face)[neighbor_index];

  return next_cell;
}

//------------------------------------------------------------------------------------------------//
/**
 * \brief Find the cell containing a position
 *
 * This function requires a linear search over the cells. You have been warned!
 *
 * \arg position. A U3 object with the appropriate dimension describing a position in
 *        space.
 *
 * \return The cell containing the point described in position, or -1 if the point is not contained
 *         by the mesh.
 */
template <int D> int CAMR_Mesh<D>::get_cell(U3 const &position) const {
  for (size_t cell = 1; cell <= active_size; ++cell) {
    double const delta = get_cell_width(cell);
    sf_double const cell_coords = get_vertex(cell);

    bool passed = true;
    for (uint32_t dimension = 1; (dimension <= D && passed); ++dimension) {
      uint32_t dimension_index = dimension - 1;
      if ((position[dimension_index] < cell_coords[dimension_index]) ||
          (position[dimension_index] >= cell_coords[dimension_index] + delta))
        passed = false;
    }
    if (passed) {
      Check(cell < INT32_MAX);
      return static_cast<int>(cell);
    }
  }
  return -1;
}

//------------------------------------------------------------------------------------------------//
/*!
 * \brief calculate the center of the portion of the face shared by two neighboring cells
 *
 * This function calculates the geometric center of the portion of the face shared by two
 * neighboring cells (i.e., the intersection of their two corresponding faces). This situations
 * arises at the boundary between refinement levels.
 *
 * \param cell cell index
 * \param face face index
 * \param next_cell neighboring cell index
 * \return face center
 */
template <int D>
U3 CAMR_Mesh<D>::get_shared_face_center(const size_t cell, const uint32_t face,
                                        const int next_cell) const {
  Require(cell > 0 && static_cast<uint32_t>(cell) <= layout->get_num_owned_cells());
  Require(face >= 1 && face <= layout->num_coarse_faces());
  Require(next_cell != 0);
  Require(next_cell != static_cast<int>(cell));

  // get the coarse face center (and just return it if the face is not refined)
  Check(cell < INT32_MAX && face < INT32_MAX);
  U3 center = get_face_center(cell, face);

  // if the face is refined, we need to calculate the center of the shared face
  if (layout->num_cells_across(static_cast<uint32_t>(cell), face) > 1) {
    // get a vector of neighboring cell indices
    sf_int const neighbors((*layout)(static_cast<uint32_t>(cell), face));

    // get the refined face index
    auto const refined_face =
        std::distance(neighbors.begin(), std::find(neighbors.begin(), neighbors.end(), next_cell)) +
        1;

    Check(static_cast<uint32_t>(refined_face) <= layout->num_amr_faces());
    Check(neighbors[refined_face - 1] == next_cell);

    // get the face vertices
    vf_double const face_vertices = get_vertices(static_cast<int>(cell), static_cast<int>(face));

    // get the dimension (x, y, or z) of the face
    Check(face < UINT32_MAX);
    int const &face_dimension = layout->face_axis_sign(static_cast<uint32_t>(face)).first;

    for (int d = 1; d <= D; ++d) {
      // for the face dimension, the coarse face center coordinate is already
      // correct
      if (d == face_dimension)
        continue;

      // for other dimensions, average the coordinates of the coarse face center and the face vertex
      // on the opposite corner of the shared face
      center[d - 1] = 0.5 * (center[d - 1] + face_vertices[d - 1][refined_face - 1]);
    }
  }

  Ensure(get_face(center, cell) == face);
  Ensure(this->next_cell(static_cast<int>(cell), static_cast<int>(face), center) == next_cell);
  return center;
}

//------------------------------------------------------------------------------------------------//
/**
 * \brief Get the distance from a point, in a given direction, to the first-crossed boundary of the
 *        given cell.
 *
 * \arg position The position
 * \arg direction The direction
 * \arg cell The cell, whose faces we are looking for crossings of.
 * \arg face. A return argument containing the face of cell which is crossed first.
 */
template <int D>
double CAMR_Mesh<D>::get_db(U3 const &position, U3 const &direction, int const cell,
                            uint32_t &face) const {
  Require(cell > 0 && static_cast<uint32_t>(cell) <= layout->get_num_owned_cells());

  sf_double cell_coords(get_vertex(cell));

  // initialize invalid face (1-based numbering)
  face = 0;
  double delta = get_cell_width(cell);
  auto winner = static_cast<double>(std::numeric_limits<float>::max());

  int positive_face;

  // Loop over the dimensions of the space
  for (uint32_t dimension = 1; dimension <= D; ++dimension) {

    uint32_t dimension_index = dimension - 1;
    const double dir = direction[dimension_index];

    // XXX would a tolerance value > 0 be better here?
    if (dir > 0.0)
      positive_face = 1;
    else if (dir < 0.0)
      positive_face = 0;
    else
      continue;

    // XXX The face 1-indexing scheme is invoked here!
    const int current_face = 2 * dimension + positive_face - 1;
    const double face_value = cell_coords[dimension_index] + positive_face * delta;
    const double distance = std::max(0.0, (face_value - position[dimension_index]) / dir);

    if (distance < winner) {
      winner = distance;
      face = current_face;
    }
  }

  Ensure(winner >= 0.0);
  Ensure(face != 0);

  return winner;
}

//------------------------------------------------------------------------------------------------//
/**
 * \brief Return the minimum distance to a cell boundary from the effective dimensions.
 *
 * \arg position A point with a cell of the mesh
 * \arg cell The cell in the mesh containing the point.
 * \return nearest orthogonal distance in the effective dimensions
 */
template <int D>
double CAMR_Mesh<D>::get_orthogonal_dist_to_bnd(const U3 &position, const int cell) const {
  Require(cell > 0 && static_cast<uint32_t>(cell) <= layout->get_num_owned_cells());

  // The maximal sphere radius is the minimum of the distances between the coordinates of the point
  // and the sides of the cell.

  sf_double cell_coords(get_vertex(cell));
  double delta = get_cell_width(cell);

  auto min_distance = static_cast<double>(std::numeric_limits<float>::max());
  for (uint32_t dimension = 1; dimension <= D; ++dimension) {

    uint32_t dimension_index = dimension - 1;
    const double low_distance = position[dimension_index] - cell_coords[dimension_index];
    const double high_distance = delta - low_distance;

    Require(low_distance >= 0.0);
    Require(high_distance >= 0.0);

    const double distance = std::min(low_distance, high_distance);
    min_distance = std::min(min_distance, distance);
  }

  return min_distance;
}

//------------------------------------------------------------------------------------------------//
/**
 * \brief Return the normal to a given face of a given cell.
 *
 * The cell argument is not used in this function, because they all have the same face normals.
 *
 * XXX This function does not use a const static array of precomputed U3's as done in
 * Cylindrical_Mesh. This would probably be a useful optimization, but initializing the arrays
 * correctly for both 2 and 3 dimensions with the available interface to U3 is not
 * obvious.
 *
 * \arg cell The cell
 * \arg face The face
 */
template <int D> U3 CAMR_Mesh<D>::get_normal(int const /*cell*/, uint32_t const face) const {
  Require(face >= 1 && face <= layout->num_coarse_faces());

  std::pair<int, bool> const face_data(layout->face_axis_sign(face));
  int const &dimension = face_data.first;
  bool const &positive = face_data.second;
  int const sign = (positive) ? 1 : -1;

  U3 normal{{0.0, 0.0, 0.0}};
  normal[dimension - 1] = sign;
  return normal;
}

//------------------------------------------------------------------------------------------------//
/**
 * \brief Return the inward normal to a given face of a given cell.
 *
 * The cell argument is not used in this function, because they all have the same face normals.
 *
 * XXX This function does not use a const static array of precomputed U3's as done in
 * Cylindrical_Mesh. This would probably be a useful optimization, but initializing the arrays
 * correctly for both 2 and 3 dimensions with the available interface to U3 is not
 * obvious.
 *
 * \arg cell The cell
 * \arg face The face
 */
template <int D> U3 CAMR_Mesh<D>::get_normal_in(size_t const /*cell*/, uint32_t const face) const {
  Require(face >= 1 && face <= layout->num_coarse_faces());

  std::pair<int, bool> const face_data(layout->face_axis_sign(face));
  int const &dimension = face_data.first;
  bool const &positive = face_data.second;
  int const sign = (positive) ? -1 : 1;

  U3 normal{{0.0, 0.0, 0.0}};
  normal[dimension - 1] = sign;
  return normal;
}

//------------------------------------------------------------------------------------------------//
/**
 * \brief Return the volume of a cell.
 * \arg cell The cell
 *
 * XXX Some way of eliminating the run-time check D==3 might help performance.
 */

// Suppress the warning: 'conditional expression is constant'. This can be removed once we adopt C++
// and the syntax is changed to 'if constexpr'
#if defined(MSVC)
#pragma warning(push)
#pragma warning(disable : 4127)
#endif

template <int D> double CAMR_Mesh<D>::volume(const size_t cell) const {
  Require(cell > 0 && static_cast<uint32_t>(cell) <= layout->get_num_owned_cells());
  double const delta = get_cell_width(cell);
  double volume = delta * delta;
  //! \bug This should be 'if constexpr' (C++17)
  if (D == 3)
    volume *= delta;
  Ensure(volume > 0);
  return volume;
}

#if defined(MSVC)
#pragma warning(pop)
#endif

//------------------------------------------------------------------------------------------------//
/**
 * \brief Return the area of a face of a cell.
 * \param[in] cell The cell.
 * \param[in] face The face.
 *
 * The face argument is not used, since all faces of a cell have the same widths.
 *
 * \note Some way of eliminating the run-time check D==3 might help performance.
 */

// Suppress the warning: 'conditional expression is constant'. This can be removed once we adopt C++
// and the syntax is changed to 'if constexpr'
#if defined(MSVC)
#pragma warning(push)
#pragma warning(disable : 4127)
#endif

template <int D> double CAMR_Mesh<D>::face_area(size_t const cell, size_t const /*face*/) const {
  Require(cell > 0 && static_cast<uint32_t>(cell) <= layout->get_num_owned_cells());
  double area = get_cell_width(cell);
  //! \bug This should be 'if constexpr'
  if (D == 3)
    area *= area;
  Ensure(area > 0);
  return area;
}

#if defined(MSVC)
#pragma warning(pop)
#endif

//------------------------------------------------------------------------------------------------//
/**
 * \brief Return the mean chord length of a cell.
 * \arg cell The cell.
 *
 * Uses Dirac's formulation l = 4*V/A for mean chord length.
 *
 * \return mean chord length
 */
template <int D> double CAMR_Mesh<D>::mean_chord_length(int const cell) const {
  Require(cell > 0 && static_cast<uint32_t>(cell) <= layout->get_num_owned_cells());

  const double mywidth = get_cell_width(cell);

  Ensure(width > 0 && D >= 1);
  return (2.0 * mywidth / static_cast<double>(D));
}

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Return the minimum distance across a cell
 *
 * \param cell CAMR_Mesh cell
 *
 * \return minimum width of the cell in any dimension
 */
template <int D> double CAMR_Mesh<D>::get_min_width(size_t const cell) const {
  Check(cell < INT32_MAX);
  return get_cell_width(static_cast<int>(cell));
}

//------------------------------------------------------------------------------------------------//
/**
 * \brief Return a vector<vector> containing the vertices of the cell.
 *
 * The data is in vertex, dimension order, i.e. inner vector is over the vertices, the outer is over
 * the dimensions.
 *
 * \arg cell The cell
 * \return A vector<vector> containing the vertices of the cell.
 */
template <int D> typename CAMR_Mesh<D>::vf_double CAMR_Mesh<D>::get_vertices(const int cell) const {

  Require(cell > 0 && static_cast<uint32_t>(cell) <= layout->get_num_owned_cells());

  const sf_double cell_vertex = get_vertex(cell);
  const int vertex_count = (D == 2 ? 4 : 8);
  const double delta = get_cell_width(cell);

  vf_double myvertices(D, sf_double(vertex_count));

  // Loop over the vertices
  for (int vertex = 0; vertex < vertex_count; ++vertex) {

    std::bitset<D> vertex_bits(vertex);

    // Loop over the dimensions.
    for (int dimension = 1; dimension <= D; ++dimension) {

      size_t dimension_index = dimension - 1;
      bool positive = vertex_bits[dimension_index];
      if (dimension_index == 0 && ((vertex > 1 && vertex < 4) || vertex > 5))
        positive = !positive;

      const double value = cell_vertex[dimension_index] + (positive ? delta : 0.0);

      myvertices[dimension_index][vertex] = value;
    }
  }
  return myvertices;
}

//------------------------------------------------------------------------------------------------//
/**
 * \brief Return a vector<vector> containing the vertices of a cell's face.
 *
 * The data is in vertex, dimension order, i.e. inner vector is over the vertices, the outer is over
 * the dimensions.
 *
 * \arg cell The cell
 * \arg face The face
 *
 * \return A vector<vector> containing the vertices of the cell.
 */

// Suppress the warning: 'conditional expression is constant'. This can be removed once we adopt C++
// and the syntax is changed to 'if constexpr'
#if defined(MSVC)
#pragma warning(push)
#pragma warning(disable : 4127)
#endif

template <int D>
typename CAMR_Mesh<D>::vf_double CAMR_Mesh<D>::get_vertices(const int cell,
                                                            const uint32_t face) const {

  Require(cell > 0 && static_cast<uint32_t>(cell) <= layout->get_num_owned_cells());
  Require(face >= 1 && face <= layout->num_coarse_faces());

  const sf_double cell_vertex = get_vertex(cell);
  const int vertex_count = (D == 2 ? 4 : 8);
  Remember(size_t const face_vertex_count = vertex_count / 2);
  const double delta = get_cell_width(cell);

  const std::pair<int, bool> face_data = layout->face_axis_sign(face);
  const int &face_dimension = face_data.first;
  const bool &face_positive = face_data.second;

  vf_double myvertices(D, sf_double());

  // Loop over the vertices
  for (int vertex = 0; vertex < vertex_count; ++vertex) {
    std::bitset<D> vertex_bits(vertex);

    // Is the vertex on the desired face?
    if (vertex_bits[face_dimension - 1] != face_positive)
      continue;

    // Loop over the dimensions.
    for (size_t dimension = 1; dimension <= D; ++dimension) {
      size_t dimension_index = dimension - 1;
      bool positive = vertex_bits[dimension_index];
      double const value = cell_vertex[dimension_index] + (positive ? delta : 0.0);
      myvertices[dimension_index].push_back(value);
    }
  }

  Ensure(myvertices[0].size() == face_vertex_count);
  Ensure(myvertices[1].size() == face_vertex_count);
  //! \bug This Ensure should be in a 'if constexpr' block.
  Ensure(D == 2 || myvertices[2].size() == face_vertex_count);
  return myvertices;
}

#if defined(MSVC)
#pragma warning(pop)
#endif

//------------------------------------------------------------------------------------------------//
/**
 * \brief Map a surface string descriptor to the numeric face value
 *
 * \arg name The string surface descriptor
 * \arg cell The cell. For interface compatibility only; Not used in this function
 * \return An integer value for the face.
 */
template <int D>
size_t CAMR_Mesh<D>::get_bndface(std::string const &name, int const /*cell*/) const {
  uint32_t face(0);
  if (name == "lox")
    face = 1;
  else if (name == "hix")
    face = 2;
  else if (name == "loy")
    face = 3;
  else if (name == "hiy")
    face = 4;
  else if (name == "loz")
    face = 5;
  else if (name == "hiz")
    face = 6;
  else
    Insist(0, "Unknown boundary descriptor string."); // LCOV_EXCL_LINE
  Require(face >= 1 && face <= layout->num_coarse_faces());
  return face;
}

//------------------------------------------------------------------------------------------------//
/**
 * \brief Return a vector of cells in this mesh which are on the given physical surface
 *
 * \arg name  The string surface descriptor
 * \return A std::vector<int> containing the cells of the mesh which have a physical boundary on
 *         the given surface.
 */
template <int D>
typename CAMR_Mesh<D>::sf_int CAMR_Mesh<D>::get_surcells(std::string const &name) const {
  Check(get_bndface(name, 1) < UINT32_MAX);
  auto const face = static_cast<uint32_t>(get_bndface(name, 1));
  Require(face >= 1 && face <= layout->num_coarse_faces());
  return layout->get_edge_cells(face);
}

//------------------------------------------------------------------------------------------------//
/**
 * \brief Verify that all of the provided cells have a vacuum boundary condition on the surface
 *        given by the string surface descriptor
 *
 * \arg name The string surface descriptor.
 * \arg cells A vector of cell indices.
 * \return boolean, true if all of the listed cells have the correct boundary condition.
 */
template <int D>
bool CAMR_Mesh<D>::check_defined_surcells(std::string const &name, sf_int const &cell_list) const {

  Check(get_bndface(name, 1) < UINT32_MAX);

  const auto face = static_cast<uint32_t>(get_bndface(name, 1));

  Require(face >= 1 && face <= layout->num_coarse_faces());

  for (auto cell : cell_list) {
    Require(cell > 0 && static_cast<uint32_t>(cell) < UINT32_MAX &&
            static_cast<uint32_t>(cell) <= layout->get_num_owned_cells());

    sf_int neighbor_data = (*layout)(static_cast<uint32_t>(cell), face);

    for (auto n : neighbor_data)
      if (n != 0)
        return false;
  }
  return true;
}

//------------------------------------------------------------------------------------------------//
/**
 * \brief Return the indices of all of the neighbors of a cell.
 * \arg cell The cell
 * \return A vector containing the indices of the neighboring cells.
 */
template <int D>
typename CAMR_Mesh<D>::sf_int CAMR_Mesh<D>::get_neighbors(uint32_t const cell) const {
  Require(cell > 0 && static_cast<uint32_t>(cell) <= layout->get_num_owned_cells());

  size_t faces = layout->num_coarse_faces();
  size_t refined_neighbors_per_face = layout->num_amr_faces();
  size_t max_neighbors = faces * refined_neighbors_per_face;

  sf_int neighbor_list;
  neighbor_list.reserve(max_neighbors);

  for (uint32_t face = 1; face <= faces; ++face) {
    sf_int neighbor_data = (*layout)(cell, face);
    std::copy(neighbor_data.begin(), neighbor_data.end(), std::back_inserter(neighbor_list));
  }

  Ensure(neighbor_list.size() <= max_neighbors);
  return neighbor_list;
}

//------------------------------------------------------------------------------------------------//
// Dump methods
//------------------------------------------------------------------------------------------------//

//------------------------------------------------------------------------------------------------//
//! Produce a vector containing a single repeated value.
template <int D> typename CAMR_Mesh<D>::sf_int CAMR_Mesh<D>::get_cell_types() const {
  //! \todo Use initializer list to avoid repeating the return type 'sf_int'
  return sf_int(active_size, (D == 2) ? rtt_viz::four_node_quadrangle : // NOLINT
                                 rtt_viz::eight_node_hexahedron);
}
//------------------------------------------------------------------------------------------------//
/**
 * \brief return a vector<vector> Victor, containing *all* of the points for *all* of the cells in
 *        the mesh.
 *
 * The data is returned in [vertex][dimension] order. That is, the dimensions of each vertex are
 * listed consecutively. Then all of the vertices of the cell are listed consecutively, within the
 * larger list of all cells.
 */
template <int D> typename CAMR_Mesh<D>::vf_double CAMR_Mesh<D>::get_point_coord() const {
  size_t const vertex_count = (D == 2 ? 4 : 8);
  size_t const mesh_size = get_size();
  size_t const result_size = mesh_size * vertex_count;

  vf_double results(result_size, sf_double(D, 0.0));
  Check(mesh_size < INT32_MAX);
  for (size_t cell = 1; cell <= mesh_size; ++cell) {
    vf_double myvertices = get_vertices(static_cast<int>(cell));
    Require(myvertices.size() == D);
    for (size_t vertex = 0; vertex < vertex_count; ++vertex) {
      size_t const result_vertex = (cell - 1) * vertex_count + vertex;
      for (size_t dim = 1; dim <= D; ++dim) {
        size_t dim_index = dim - 1;
        results[result_vertex][dim_index] = myvertices[dim_index][vertex];
      }
    }
  }
  return results;
}

//------------------------------------------------------------------------------------------------//
/**
 * \brief Return an index map from cells and vertices to the entries in the table returned from
 *        get_point_coord().
 *
 * e.g.
 * \code
 *    vector<vector<double> > vertices = mesh.get_point_coord()
 *    vector<vector<int> >    indices  = mesh.get_cell_pair()
 * \endcode
 *
 * Then indices[cell-1][vertex] gives the first index into vertices for the coordinates of that
 * vertex.
 */
template <int D> typename CAMR_Mesh<D>::vf_int CAMR_Mesh<D>::get_cell_pair() const {
  int const vertex_count = (D == 2 ? 4 : 8);
  size_t const mesh_size = get_size();

  vf_int indices(mesh_size, sf_int(vertex_count));

  size_t index = 0;
  for (size_t cell = 1; cell <= mesh_size; ++cell) {
    for (size_t vertex = 0; vertex < vertex_count; ++vertex) {
      indices[cell - 1][vertex] = static_cast<int>(++index);
    }
    Require(index % vertex_count == 0);
    Require(index / vertex_count == cell);
  }

  Ensure(index == mesh_size * vertex_count);
  return indices;
}

//------------------------------------------------------------------------------------------------//
//! Comparison operator
template <int D> bool CAMR_Mesh<D>::operator==(CAMR_Mesh<D> const &rhs) const {
  return (layout == rhs.layout) && (vertices == rhs.vertices) && (levels == rhs.levels);
}

} // namespace rtt_mesh

#endif // rtt_mesh_CAMR_Mesh_t_hh

//------------------------------------------------------------------------------------------------//
// end of mesh/CAMR_Mesh.t.hh
//------------------------------------------------------------------------------------------------//
