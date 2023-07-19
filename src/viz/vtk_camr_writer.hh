//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   viz/vtk_camr_writer.hh
 * \author H. Park
 * \date   Mon Jun 19 2023
 * \brief  vtk writer for CAMR file
 * \note   Copyright (C) 2023 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#ifndef rtt_viz_vtk_camr_writer_hh
#define rtt_viz_vtk_camr_writer_hh

#include <sstream>
#include <type_traits>
#include <unordered_map>
#include <vector>

namespace rtt_viz {

//================================================================================================//
/*!
 * \class vtk_camr_writer_t
 *
 * \brief A VTK writer class for legacy files.
 */
//================================================================================================//
class vtk_camr_writer_t {

public:
  /// constructors
  vtk_camr_writer_t() = default;
  explicit vtk_camr_writer_t(int /*rank*/) {}

  /// Open a file
  int open(const char *filename) {
    file_.open(filename);
    file_ << "<?xml version=\"1.0\"?>" << std::endl;

    file_.setf(std::ios::scientific);
    file_.precision(sigfigs_);

    return !file_.good();
  }

  /// Close a file
  int close() {
    file_.close();
    return !file_.good();
  }

public:
  // write "VTK files "
  int write_camr_unstructured_begin() {
    file_ << R"(<VTKFile type="UnstructuredGrid" version="2.0" byte_order="LittleEndian">)"
          << std::endl;
    file_ << "<UnstructuredGrid>" << std::endl;

    return !file_.good();
  }
  int write_camr_unstructured_end() {
    file_ << "</UnstructuredGrid>" << std::endl;
    file_ << "</VTKFile>" << std::endl;
    return !file_.good();
  }

  //------------------------------------------------------------------------------------------------//
  // Write piece data
  int write_camr_piece_begin(const int num_dims, const int ncells) {
    int nvts = static_cast<int>(pow(2, num_dims)) * ncells;
    file_ << R"(<Piece NumberOfPoints=")" << nvts << R"(" NumberOfCells=")" << ncells << R"(">)"
          << std::endl;
    return !file_.good();
  }

  int write_camr_piece_end() {
    file_ << "</Piece>" << std::endl;
    return !file_.good();
  }

  //------------------------------------------------------------------------------------------------//
  // Write point coordinates
  int write_camr_points(const std::size_t num_dims, const std::size_t ncells, const double dx,
                        const std::vector<double> &vertex, const std::vector<int> &levels) {

    file_ << "<Points>" << std::endl;
    file_ << R"(<DataArray type="Float64" NumberOfComponents="3">)" << std::endl;
    ;

    for (int i = 0; i < static_cast<int>(ncells); ++i) {
      double width = dx / pow(2, levels[i]);
      double x0 = vertex[i * num_dims];
      double y0 = vertex[i * num_dims + 1];
      double x1 = x0 + width;
      double y1 = y0 + width;
      if (num_dims == 3) {
        double z0 = vertex[i * num_dims + 2];
        double z1 = z0 + width;

        file_ << x0 << " " << y0 << " " << z0 << std::endl;
        file_ << x1 << " " << y0 << " " << z0 << std::endl;
        file_ << x1 << " " << y1 << " " << z0 << std::endl;
        file_ << x0 << " " << y1 << " " << z0 << std::endl;
        file_ << x0 << " " << y0 << " " << z1 << std::endl;
        file_ << x1 << " " << y0 << " " << z1 << std::endl;
        file_ << x1 << " " << y1 << " " << z1 << std::endl;
        file_ << x0 << " " << y1 << " " << z1 << std::endl;

      } else {
        file_ << x0 << " " << y0 << " " << 0.0 << std::endl;
        file_ << x1 << " " << y0 << " " << 0.0 << std::endl;
        file_ << x1 << " " << y1 << " " << 0.0 << std::endl;
        file_ << x0 << " " << y1 << " " << 0.0 << std::endl;
      }

    } //i
    file_ << "</DataArray>" << std::endl;

    file_ << "</Points>" << std::endl;

    return !file_.good();
  }

  //------------------------------------------------------------------------------------------------//
  // Write cells either quad (2d) or hex (3d)
  int write_camr_cells(const int num_dims, const int ncells) {
    int nvts = static_cast<int>(pow(2, num_dims));
    int element_type = num_dims == 2 ? 9 : 12;

    file_ << "<Cells>" << std::endl;
    // connectivity.
    file_ << R"(<DataArray type="Int32" Name="connectivity">)" << std::endl;
    ;
    int cntr(0);

    for (int i = 0; i < ncells; ++i) {

      for (int j = 0; j < nvts; ++j)
        file_ << cntr++ << " ";

      file_ << std::endl;
    }
    file_ << "</DataArray>" << std::endl;

    // offsets
    cntr = nvts;
    file_ << R"(<DataArray type="Int32" Name="offsets">)" << std::endl;
    ;
    for (int i = 0; i < ncells; ++i) {
      file_ << cntr << " ";
      cntr += nvts;
    }
    file_ << std::endl;
    file_ << "</DataArray>" << std::endl;

    // types
    file_ << R"(<DataArray type="UInt8" Name="types">)" << std::endl;
    ;
    for (int i = 0; i < ncells; ++i) {
      file_ << element_type << " ";
    }
    file_ << std::endl;

    file_ << "</DataArray>" << std::endl;

    file_ << "</Cells>" << std::endl;

    return !file_.good();
  }

  //------------------------------------------------------------------------------------------------//
  // Write cell data
  int write_camr_cell_data_begin() {
    file_ << "<CellData>" << std::endl;
    return !file_.good();
  }

  int write_camr_cell_data_end() {
    file_ << "</CellData>" << std::endl;
    return !file_.good();
  }

  template <typename T>
  int write_camr_cell_data(const std::string &name, const T begin, const T end) {
    file_ << R"(<DataArray type="Float64" Name=")" << name << R"(")"
          << R"( NumberOfComponents="1">)" << std::endl;

    for (auto i = begin; i != end; ++i)
      file_ << *i << " ";

    file_ << "\n</DataArray>\n";

    return !file_.good();
  }

  int write_camr_cell_data(const std::string &name, const std::vector<double> &data) {
    int data_size = static_cast<int>(data.size());
    file_ << R"(<DataArray type="Float64" Name=")" << name << R"(">)" << std::endl;
    for (int i = 0; i < data_size; ++i)
      file_ << data[i] << " ";
    file_ << std::endl;
    file_ << "</DataArray>" << std::endl;
    return !file_.good();
  }

  //------------------------------------------------------------------------------------------------//
  // Write root file
  int write_camr_pvtu_file(const int nranks, const int /*num_dims*/, const double time,
                           const std::string &file_prefix,
                           const std::vector<std::string> &var_names) {
    file_ << R"(<VTKFile type="PUnstructuredGrid" version="2.0" byte_order="LittleEndian">)"
          << std::endl;
    file_ << R"(<PUnstructuredGrid GhostLevel="0">)" << std::endl;
    file_ << "<FieldData>" << std::endl;
    file_ << R"(<DataArray type="Float64" Name="TimeValue" NumberOfTuples="1">)" << time
          << std::endl;
    file_ << "</DataArray>" << std::endl;
    file_ << "</FieldData>" << std::endl;

    for (int i = 0; i < nranks; ++i)
      file_ << R"(<Piece Source=")" << file_prefix << i << R"().vtu"/>)" << std::endl;

    file_ << "<PPoints>" << std::endl;
    file_ << R"(<PDataArray type="Float64" NumberOfComponents="3"/>)" << std::endl;
    file_ << "</PPoints>" << std::endl;
    file_ << "<PCells>" << std::endl;

    file_ << R"(<PDataArray type="Int32" Name="connectivity"/>)" << std::endl;
    file_ << R"(<PDataArray type="Int32" Name="offsets"/>)" << std::endl;
    file_ << R"(<PDataArray type="UInt8" Name="types"/>)" << std::endl;

    file_ << "</PCells>" << std::endl;

    file_ << "<PCellData>" << std::endl;

    for (auto &var_name : var_names)
      file_ << R"(<PDataArray type="Float64"  Name=")" << var_name
            << R"(" NumberOfComponents="1"/>)" << std::endl;
    file_ << "</PCellData>" << std::endl;

    file_ << "</PUnstructuredGrid>" << std::endl;

    file_ << "</VTKFile>" << std::endl;

    return !file_.good();
  }

private:
  //! file pointer
  std::ofstream file_;

  //! significant figures
  int sigfigs_ = std::numeric_limits<double>::digits10;
};

} // namespace rtt_viz

#endif // rtt_viz_vtk_camr_writer_hh

//------------------------------------------------------------------------------------------------//
// end of viz/vtk_camr_writer.hh
//------------------------------------------------------------------------------------------------//
