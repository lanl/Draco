//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   mc/OS_Builder.hh
 * \author Thomas M. Evans
 * \date   Mon Feb  9 16:16:07 1998
 * \brief  OS_Builder class header file.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef __mc_OS_Builder_hh__
#define __mc_OS_Builder_hh__

#include "Coord_sys.hh"
#include "Layout.hh"
#include "OS_Mesh.hh"
#include "ds++/SP.hh"
#include "ds++/Assert.hh"
#include <vector>
#include <string>
#include <iostream>
#include <fstream>

namespace rtt_mc 
{

//===========================================================================//
/*!
 * \class OS_Builder
 *
 * The OS_Builder class builds an instance of rtt_mc::OS_Mesh using a simple
 * MC-defined orthogonal mesh format.  Future incantations may use more
 * advanced mesh format readers in order to build OS meshes from more
 * advanced formats.
 * 
 * /sa The examples page (mc/test/tstOSMesh.cc) for examples how to build OS
 * meshes using OS_Builder. 
 */
// revision history:
// -----------------
//  0) original
//  1)  3-18-98 : added generalized mesh constructor which consists of 
//                calculating vertex-based arrays and sending them to the 
//                OS_Mesh constructor
//  2)   4-6-99 : made OS_Builder class templated on an interface type
//  3)  4-13-99 : moved into mc package
//  4) 10-FEB-00: OS_Builder now reads its own mesh format file.  All it
//                requires to do this is the name of the input file.  The
//                input file can be the same as the regular Milagro input
//                file
//===========================================================================//

class OS_Builder
{
  public:
    // Useful typdefs to std:: namespace members.
    typedef rtt_dsxx::SP<OS_Mesh>             SP_Mesh;
    typedef rtt_dsxx::SP<Coord_sys>           SP_Coord_sys;
    typedef rtt_dsxx::SP<Layout>              SP_Layout;
    typedef std::vector<int>                  sf_int;
    typedef std::vector<std::vector<int> >    vf_int;
    typedef std::vector<double>               sf_double;
    typedef std::vector<std::vector<double> > vf_double;
    typedef std::vector<std::string>          sf_string;
    typedef std::string                       std_string;
    typedef std::ifstream                     std_ifstream;
    
  private:
    // Data from Parser needed to build mesh.
    
    // Input file.
    std_string mesh_file;

    // Coordinate system string.
    std_string coord_system;

    // Number of fine cells per coarse cell.
    vf_int fine_cells;
    
    // Ratio for fine zoning within a coarse cell (next fine cell 
    // is "ratio" times larger/smaller.)  Either all default to 
    // unity or all must be entered.
    vf_double fine_ratio;

    // Recursive total number of fine_cells per coarse cell.
    vf_int accum_cells;

    // Coarse edges.
    vf_double coarse_edge;

    // Number of fine_cells along each dimension.
    vf_double fine_edge;

    // Boundary conditions.
    sf_string bnd_cond;

    // Zone map.
    sf_int zone;
    vf_int cell_zone;
  
    // Defined cell regions.
    vf_int regions;

    // Surface source positional information.
    sf_string ss_pos;
    sf_int    num_defined_surcells;
    vf_int    defined_surcells;

    // Pointer to built Mesh.
    SP_Mesh mesh;

    // Member functions for building OS_Mesh

    // Parse the mesh input file.
    void parser();
    void parser2D(std_ifstream &);
    void parser3D(std_ifstream &);
    void source_parser(std_ifstream &);

    // Build Layout helper functions.
    SP_Layout build_Layout(const Coord_sys &);
    void assign2D(Layout &);
    void assign3D(Layout &);

    // Build Coord_sys helper functions.
    SP_Coord_sys build_Coord();

    // Build Mesh helper functions.
    SP_Mesh build_2DMesh(SP_Coord_sys, Layout &);
    SP_Mesh build_3DMesh(SP_Coord_sys, Layout &);

    // Member functions for cell-zone mapping
    void zone_mapper();
    void cell_zoner(int, int);
    void cell_zoner(int, int, int);

    // Calculate defined surface cells.
    void calc_defined_surcells();

  public:
    // Constructor.
    template<class IT> explicit OS_Builder(rtt_dsxx::SP<IT>);

    // Build Mesh function.
    SP_Mesh build_Mesh();

    // Map a zone centered field into a cell centered field.
    template<class T> 
    std::vector<T> zone_cell_mapper(const std::vector<T> &) const;

    // ACCESSORS
    
    // Get a copy of the built mesh.
    SP_Mesh get_Mesh() const { Require(mesh); return mesh; }

    // Get cell regions for graphics dumping.
    sf_int get_regions() const;
    int get_num_regions() const { return regions.size(); }

    // Get cell zone information.
    int get_num_zones() const { return cell_zone.size(); }
    sf_int get_cells_in_zone(int z) const { return cell_zone[z-1]; }

    // Get the defined surcells list and positions.
    vf_int get_defined_surcells() const;
    sf_string get_ss_pos() const { return ss_pos; }
};

//---------------------------------------------------------------------------//
// Templated functions for OS_Builder
//---------------------------------------------------------------------------//
// Constructor.

template<class IT>
OS_Builder::OS_Builder(rtt_dsxx::SP<IT> interface)
    : mesh_file(),
      coord_system(), 
      fine_cells(),
      fine_ratio(),
      accum_cells(),
      coarse_edge(),
      fine_edge(),
      bnd_cond(),
      zone(),
      cell_zone(),
      regions(),
      ss_pos(),
      num_defined_surcells(),
      defined_surcells(),
      mesh()
{
    Require (interface);

    // get mesh input file name from interface
    mesh_file = interface->get_mesh_file();

    // parse the mesh input file
    parser();

    // do zone mapping
    zone_mapper();

    Ensure (!mesh);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Map zone centered data fields into cell centered data fields.
 *
 * This function takes a zone centered field of size [0:Nz-1] where Nz is the
 * number of zones in the problem, and it maps the data into a cell centered
 * field of size [0:Nc-1].  Here Nc is the number of mesh cells in the
 * problem.  The builder has intrinsically (after construction) the
 * cell-to-zone mappings that make this operation possible.  The field types
 * must be vectors of type T.
 *
 * \param zone_field zone centered field to be converted into cell centered
 * field. 
 */
template<class T>
std::vector<T> OS_Builder::zone_cell_mapper(const std::vector<T> &zone_field)
    const 
{
    // we will use vector throughout this function
    using std::vector;

    Require (zone_field.size() == cell_zone.size());

    // make a cell-sized return vector
    vector<T> cell_field(zone.size());

    // assign cell values to cell_field based on zonal values
    for (int cell = 0; cell < cell_field.size(); cell++)
	cell_field[cell] = zone_field[zone[cell]-1];

    // return the mapped field
    return cell_field;
}

} // end namespace rtt_mc

#endif                          // __mc_OS_Builder_hh__

//---------------------------------------------------------------------------//
//                              end of mc/OS_Builder.hh
//---------------------------------------------------------------------------//
