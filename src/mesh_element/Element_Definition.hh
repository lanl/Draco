//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   mesh_element/Element_Definition.hh
 * \author John McGhee
 * \date   Fri Feb 25 10:03:18 2000
 * \brief  Header file for the RTT Element_Definition class.
 * \note   Copyright (C) 2010-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#ifndef rtt_mesh_element_Element_Definition_hh
#define rtt_mesh_element_Element_Definition_hh

#include "ds++/Assert.hh"
#include <ostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace rtt_mesh_element {

//================================================================================================//
/*!
 * \class Element_Definition
 *
 * \brief Provides some descriptive information on the standard mesh elements used in the RTT
 *        meshReader class.
 *
 * A few high points, trying not to wax eloquent. It was originally desired to create a simple class
 * that would concisely, unambiguously, and completely describe any mesh element that could be
 * conceived. While this may be a laudable goal, it appears to be harder than it appears. Perhaps we
 * could get some help on this from some computational geometry experts at some time. In the mean
 * time here is my 80% solution.
 *
 * First, we will reduce the scope from any element to just the elements currently supported by the
 * <a href="http://www.cgns.org/"> CGNS </a> data storage system. CGNS is an emerging industry
 * standard for the storage and retrival of computational physics data. Remember that it is only
 * necessary to describe the problem "geometry" with these elements, not the solution, or any other
 * field on the mesh, so this may not be as much of a restriction as it first appears. The CGNS set
 * consists of 18 elements including most of the commonly used ones. Moreover, remember we currently
 * have no means of generating a mesh with weird custom elements. Any mesh anyone in the group has
 * ever run on can be expressed with just six of the 18 CGNS elements.
 *
 * Second, we will not try to design a completely general element description, but will settle for
 * providing a limited set of services that can be used to discover a lot of things about the
 * elements in the CGNS sub-set, but may not necessarily be a universal, complete, and unambiguous
 * description. The ultimate authority on the element descriptions are the <a
 * href="http://www.CGNS.org/documents/Elements.pdf"> figures </a> and text found in the CGNS
 * SIDS-Additions manual.
 *
 * The description implemented herein utilizes a hierarchical approach. 3D elements are described as
 * assemblies of 2D elements, which are composed of 1D elements, which are themselves composed of
 * nodes. For example, a 3D hexahedra is described in terms of its 2D quadrilateral faces, which are
 * described in terms of 1D line edge elements, which are then described in terms of their
 * constituent nodes.  This approach appears to be adequate for the subset of elements under
 * consideration herein, but it is not clear that this will suffice in the general case.
 *
 * Utilities are provided to inquired about the type of a face (i.e. quad, triangle, etc....) as
 * well as the nodes that compose the face.
 *
 * In addition to face types, there is a concept of "node-location" within the element. All nodes
 * are given a location (i.e. "CORNER", "EDGE", etc....) to aide in the description of the
 * element. Again this appears to be adequate for the sub-set of elements under consideration herein
 * but may not be adequate in a more general case.
 *
 * It is hoped that the node-location, together with the data available through recursively
 * descending through element faces, edges, and nodes provides an adequate amount of information for
 * our present needs. However, it is difficult to show that this description is complete and
 * unambiguous.
 *
 * \sa The \ref overview_meshreaders page provides an overview of the other utilities in the
 *     rtt_mesh_element namespace.
 */
//================================================================================================//

class Element_Definition {

  // NESTED CLASSES AND TYPEDEFS

public:
  /*!
   * \brief Describes the location of a node within an element.
   *
   * For the purposes of this enumeration, the terms have the following meaning: A "corner" node
   * terminates one or more edges of an element.  The term "edge" describes a node that lies on the
   * interior of a 1D element. The term "face" describes a node that lies on the interior of a 2D
   * element. Finally the term "cell" connotates a node that lies on the interior of a 3D element.
   *
   * All elements will always have corner nodes. In addition, all elements may have edge nodes. Two
   * and three dimensional elements may also have face nodes, and finally, three-dimensional
   * elements may have cell nodes. Under these definitions, note that a node's location is unchanged
   * in an element and all its sub-elements. i.e. the corner nodes of a quadrilateral are also
   * corner nodes in the line elements which form the edges of the quadrilateral.
   */
  enum Node_Location {
    CORNER, /*!< Terminates one or more edges of an element. */
    EDGE,   /*!< Located in the interior of a 1D element. */
    FACE,   /*!< Located in the interior of a 2D element. */
    CELL    /*!< Located in the interior of a 3D element. */
  };

  /*!
   * \brief Standard element identifiers.
   *
   * These names and the elements that they represent are the same as those defined in the <a
   * href="http://www.cgns.org/"> CGNS </a> SIDS Manual.  <a
   * href="http://www.CGNS.org/documents/Elements.pdf"> Element-Descriptions </a> (Adobe PDF format)
   * are are available at the CGNS www site.
   */
  enum Element_Type {
    NODE,       //!< A dimensionless point in space.
    BAR_2,      //!< The basic one-D, two-node "line" element.
    BAR_3,      //!< Same as "BAR_2" except that a node is added in the center.
    TRI_3,      //!< The basic two-D, three-node, "triangle" element.
    TRI_6,      /*!< Same as "TRI_3" except that nodes are added in the * middle of each edge. This
                 *   is the standard quadratic-serendipity finite element triangle. */
    QUAD_4,     //!< The basic two-D, four-node "quadrilateral" element.
    QUAD_5,     //!< A quad with a node in the center of one face.
    QUAD_6,     /*!< A quad with nodes in the center of two ADJOINING faces. This is the default
                 *   QUAD_6. */
    QUAD_6a,    //!< A quad with nodes in the center of two ADJOINING faces.
    QUAD_6o,    //!< A quad with nodes in the center of two OPPOSITE faces.
    QUAD_7,     //!< A quad with nodes in the center of three faces.
    QUAD_8,     /*!< A quad with nodes in the center of all four faces. This is standard
                 *   quadratic-serendipity finite element quad.*/
    QUAD_9,     //!< Same as "QUAD_8" except a node is added in the center of the quad.
    PENTAGON_5, /*!< The basic two-D, five-node "pentagon" element.  Elements with this topology are
                 *   quite common in an AMR mesh. */
    HEXAGON_6,  /*!< The basic two-D, six-node "hexagon" element.  Elements with this topology are
                 *   quite common in an AMR mesh. */
    HEPTAGON_7, /*!< The basic two-D, seven-node "heptagon" element.  Elements with this topology
                 *   can occur in an AMR mesh. */
    OCTAGON_8,  /*!< The basic two-D, eight-node "octagon" element.  Elements with this topology can
                 *   occur in an AMR mesh. */
    TETRA_4,    //!< The basic three-D, four-node "tetrahedral" element.

    TETRA_10,   /*!< Same as "TETRA_4" except that a node is added in the middle of each edge. This
                 *   is the standard quadratic-serendipity finite element tet.*/
    PYRA_5,     /*!< The basic three-D, five-node, "pyramid" element.  This is a hex with one face
                 *   collapsed to a point.*/
    PYRA_14,    /*!< Same as "PYRA_5" except that a node is added on each edge, and one at the
                 *   center. */
    PENTA_6,    /*!< The basic three-D, six-node "pentahedron". Also known as a "triangular-prism",
                 *   or "wedge". */
    PENTA_15,   /*!< Same as "PENTA-6" except that nodes are added in the center of each edge. This
                 *   is the standard quadratic-serendipity finite element wedge.*/
    PENTA_18,   /*!< Same as "PENTA-15" except that nodes are added in the center of each
                 *   quadrilateral face. */
    HEXA_8,     //!< The basic three-D, eight-node "hexahedron".
    HEXA_20,    /*!< Same as "HEXA_8" except that a node is added in the center of each edge. This
                 *   is the standard quadratic-serendipity finite element hex.*/
    HEXA_27,    /*!< Same as "HEXA_20" except that a node is added in the center of each face, and
                 *   at the center of the element. */
    POLYHEDRON, //!< A hexahedron with, possibly, subdivided hexadedral neighbors.
    POLYGON,    //!< A polygon element with straight sides.
    NUMBER_OF_ELEMENT_TYPES
  };

private:
  // DATA

  std::string name;
  Element_Type type;
  unsigned dimension;
  unsigned number_of_nodes;
  unsigned number_of_sides;
  std::vector<Element_Definition> elem_defs;
  std::vector<unsigned> side_type;
  std::vector<std::vector<unsigned>> side_nodes;

public:
  // CREATORS

  /*!
   * \brief Constructor for the Element_Definition class.
   * \param type_ The element type to be constructed.
   */
  explicit Element_Definition(Element_Type const &type_);

  /*!
   * \brief Constructor for the Element_Definition class.
   *
   * This constructor supports the description of a nonstandard element type.
   *
   * \param name_ The name of the element.
   *
   * \param dimension_ The dimension of the element. i.e. nodes return 0, lines return 1, quads
   *        return 2, hexahedra return 3.
   *
   * \param number_of_nodes_ Total number of nodes in the element
   *
   * \param number_of_sides_ The number of n-1 dimensional entities that compose an n dimensional
   *        element. i.e. nodes return 0, lines return 2, quads return 4, hexahedra return 6.
   *
   * \param elem_defs_ Element definitions that describe element sides.  There need be only one such
   *        definition for each type of side present in the element.  For example, a QUAD_4 element
   *        would need only one side element definition, for BAR_2.
   *
   * \param side_type_ Index into \c elem_defs_ of the element definition appropriate for each side.
   *
   * \param side_nodes_ A vector of vectors specifying the nodes associated with each side. For
   *        example, <code>side_nodes_[2]</code> is a vector specifying the nodes associated with
   *        the third side of the element.  Note that the node numbering is 0 based.
   *
   * \pre <code>dimension_>=0</code>
   *
   * \pre <code>number_of_nodes_>0</code>
   *
   * \pre <code>number_of_sides_>=0</code>
   *
   * \pre All elements of \c elem_defs_ must satisfy
   *        <code>elem_defs_[i].get_dimension()+1==dimension_</code>
   *
   * \pre <code>side_type_.size()==number_of_sides_</code>
   *
   * \pre All elements of \c side_type_ must satisfy
   *        <code>static_cast<unsigned>(side_type_[i])<elem_defs_.size()</code>
   *
   * \pre <code>side_nodes_.size()==number_of_sides_</code>
   *
   * \pre All elements of \c side_nodes_ must satisfy
   *        <code>side_nodes_[i].size() == elem_defs_[side_type_[i]].get_number_of_nodes() </code>
   *
   * \pre All elements of \c side_nodes_ must satisfy
   *        <code>static_cast<unsigned>(side_nodes_[i][j])<number_of_nodes_</code>
   *
   * \post <code> get_type()==Element_Definition::POLYGON </code>
   *
   * \post <code> get_name()==name_  </code>
   *
   * \post <code> get_dimension()==dimension_  </code>
   *
   * \post <code> get_number_of_nodes()==number_of_nodes_  </code>
   *
   * \post <code> get_number_of_sides()==number_of_sides_  </code>
   *
   * \post <code> get_side_type(i)==elem_defs_[side_type_[i]]  </code>
   *
   * \post <code> get_side_nodes(i)==side_nodes_[i]  </code>
   */
  Element_Definition(std::string const &name_, unsigned dimension_, unsigned number_of_nodes_,
                     unsigned number_of_sides_, std::vector<Element_Definition> const &elem_defs_,
                     std::vector<unsigned> const &side_type_,
                     std::vector<std::vector<unsigned>> const &side_nodes_);

  // MANIPULATORS

  /*!
   * \brief Destructor for the Element_Definition class.
   *
   * This destructor is virtual, implying that Element_Definition is extensible by inheritance.
   */
  virtual ~Element_Definition() = default;

  //! Copy constructor
  Element_Definition(Element_Definition const &rhs) = default;
  Element_Definition(Element_Definition &&rhs) = default;

  //! Assignment operator
  Element_Definition &operator=(Element_Definition const &rhs) = default;
  Element_Definition &operator=(Element_Definition &&rhs) = default;

  // ACCESSORS

  /*!
   * \brief Returns the name of an element.
   * \return Returns the element name as a string.
   */
  std::string get_name() const { return name; }

  /*!
   * \brief Returns the type of an element.
   * \return Returns the element type.
   */
  Element_Type get_type() const { return type; }

  /*!
   * \brief Returns the total number of nodes in an element.
   * \return Total number of nodes in an element.
   */
  unsigned get_number_of_nodes() const { return number_of_nodes; }
  /*!
   * \brief Returns the dimension of an element. i.e. nodes return 0, lines return 1, quads return
   *        2, hexahedra return 3.
   *
   * \return The element dimension (0, 1, 2, or 3).
   */
  unsigned get_dimension() const { return dimension; }
  /*!
   * \brief Returns the number of sides on an element.
   *
   * \return The number of n-1 dimensional entities that compose an n dimensional
   *         element. i.e. nodes return 0, lines return 2, quads return 4, hexahedra return 6.
   */
  unsigned get_number_of_sides() const { return number_of_sides; }

  /*!
   * \brief Returns the type (i.e. quad, tri, etc.) of a specified element side.
   *
   * \return Returns a valid element definition that describes a element side. Can be queried using
   *         any of the accessors provided in the Element_Definition class.

   * \param side_number Side number for which a type is desired.  Side numbers
   *        are in the range [0:number_of_sides).
   *
   * Note that there is no valid side number for a "NODE" element.  "Side" in the context of this
   * method means the (n-1) dimensional element that composes a n dimensional element.
   */
  Element_Definition get_side_type(unsigned const side_number) const {
    Insist(side_number < side_type.size(), "Side index out of range!");
    return elem_defs[side_type[side_number]];
  }

  /*!
   * \brief Returns a vector of node numbers that are associated with a particular element side.
   *
   * \param side_number The number of the element side for which the nodes are desired. Side numbers
   *        are in the range [0:number_of_sides)
   *
   * \return A vector of the nodes associated with the side. Note that the node numbering is 0
   *         based.
   *
   * "Side" in the context of this method means the (n-1) dimensional element that composes a (n)
   * dimensional element. For example, on a hexahedra, a side is a quadrilateral, whereas, on a
   * quadrilateral a side is a line element.  The returned order of the side nodes is
   * significant. The side-node numbers are returned in the following order based on node location:
   * (corners, edges, faces, cells). For sides which are faces of 3D elements, the vector cross
   * product of the vector from (side-node1 to side-node2) with the vector from (side-node1 to side-
   * node3) results in a vector that is oriented outward from the parent element.  Equivalently, the
   * side corner-nodes are listed sequentially in a counter-clockwise direction when viewed from
   * outside the element. Both corner and edge nodes are returned in a sequential order as one
   * progresses around a side. Moreover, the corner and edge nodes are returned so that edge-node1
   * lies between corner-node1 and corner-node2, etc., etc.
   *
   * For sides which are edges of 2D elements, the vector cross product of the vector from
   * (side-node1 to side-node2) with a vector pointing towards the observer results in a vector that
   * is oriented outward from the parent element.
   *
   * Note that there is no valid side number for a "NODE" element.
   */
  std::vector<unsigned> get_side_nodes(unsigned const side_number) const {
    Insist(side_number < side_nodes.size(), "Side index out of range!");
    return side_nodes[side_number];
  }

  std::vector<unsigned> get_number_of_face_nodes() const {
    std::vector<unsigned> number_of_face_nodes(side_type.size());

    for (unsigned s = 0; s < side_type.size(); ++s)
      number_of_face_nodes[s] = get_side_type(s).get_number_of_nodes();

    return number_of_face_nodes;
  }

  std::vector<std::vector<unsigned>> get_face_nodes() const {
    std::vector<std::vector<unsigned>> face_nodes(side_nodes.size());

    for (unsigned s = 0; s < face_nodes.size(); ++s) {
      std::vector<unsigned> nodes(get_side_nodes(s));
      face_nodes[s].resize(nodes.size());

      for (unsigned n = 0; n < nodes.size(); ++n)
        face_nodes[s][n] = nodes[n];
    }
    return face_nodes;
  }

  /*!
   * \brief Performs some simple sanity checks on the private data of the Element_Description
   *        class. Note that this only works with DBC turned on.
   */
  bool invariant_satisfied() const;

  //! Prints the element description.
  std::ostream &print(std::ostream &os_out) const;

  //! Define convenience ostream inserter.
  friend std::ostream &operator<<(std::ostream &os, Element_Definition const &rhs) {
    return rhs.print(os);
  }

private:
  // IMPLEMENTATION

  void construct_node();
  void construct_bar();
  void construct_tri();
  void construct_quad();
  // void construct_pentagon();
  void construct_tetra();
  void construct_pyra();
  void construct_penta();
  void construct_hexa();
};

} // end namespace rtt_mesh_element

#endif // rtt_mesh_element_Element_Definition_hh

//------------------------------------------------------------------------------------------------//
// end of mesh_element/Element_Definition.hh
//------------------------------------------------------------------------------------------------//
