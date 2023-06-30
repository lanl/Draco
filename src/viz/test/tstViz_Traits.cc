//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   viz/test/tstViz_Traits.cc
 * \author Thomas M. Evans
 * \date   Fri Jan 21 17:51:52 2000
 * \brief  Viz_Traits test.
 * \note   Copyright (C) 2014-2023 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "ds++/Release.hh"
#include "ds++/ScalarUnitTest.hh"
#include "ds++/Soft_Equivalence.hh"
#include "viz/Viz_Traits.hh"

using namespace std;
using rtt_viz::Viz_Traits;

//------------------------------------------------------------------------------------------------//
// simple test field class for checking viz traits

template <typename T> class Test_Field {
public:
  using value_type = vector<T>;

private:
  vector<vector<T>> data;

public:
  explicit Test_Field(vector<vector<T>> data_in) : data(std::move(data_in)) {}

  T operator()(size_t i, size_t j) const { return data[i][j]; }
  size_t nrows() const { return data.size(); }
  size_t ncols(size_t r) const { return data[r].size(); }
  size_t size() const { return data.size(); }
  vector<T> operator[](size_t i) const { return data[i]; }
};

//------------------------------------------------------------------------------------------------//
// Use soft_equiv for floating-point types, but not for integral types.

bool compare_vdf_field(double const &v1, double const &v2) { return rtt_dsxx::soft_equiv(v1, v2); }
bool compare_vdf_field(float const &v1, float const &v2) { return rtt_dsxx::soft_equiv(v1, v2); }
bool compare_vdf_field(int const &v1, int const &v2) { return v1 == v2; }

//------------------------------------------------------------------------------------------------//
// test vector traits specialization

template <typename VVF> void test_vector(rtt_dsxx::UnitTest &ut) {
  using VVFet = typename VVF::value_type::value_type;

  VVF field(3);

  for (size_t i = 0; i < field.size(); i++) {
    field[i].resize(i + 2);
    for (size_t j = 0; j < field[i].size(); j++)
      field[i][j] = static_cast<VVFet>(2 * i + 4 * j);
  }

  Viz_Traits<VVF> vdf(field);

  FAIL_IF_NOT(vdf.nrows() == field.size());
  for (size_t i = 0; i < vdf.nrows(); i++) {
    FAIL_IF_NOT(vdf.ncols(i) == field[i].size());
    for (size_t j = 0; j < vdf.ncols(i); j++) {
      FAIL_IF_NOT(compare_vdf_field(vdf(i, j), field[i][j]));
      FAIL_IF_NOT(compare_vdf_field(vdf(i, j), static_cast<VVFet>(2 * i + 4 * j)));
    }
  }
  if (ut.numFails == 0)
    PASSMSG("test_vector passes.");
  return;
}

//------------------------------------------------------------------------------------------------//
// standard Viz_Traits field test

template <typename T> void test_FT(rtt_dsxx::UnitTest &ut) {
  vector<vector<T>> field(3);
  for (size_t i = 0; i < field.size(); i++) {
    field[i].resize(i + 2);
    for (size_t j = 0; j < field[i].size(); j++)
      field[i][j] = static_cast<T>(2 * i + 4 * j);
  }

  Test_Field<T> test_field(field);
  Viz_Traits<Test_Field<T>> vt(test_field);

  FAIL_IF_NOT(vt.nrows() == 3);
  for (size_t i = 0; i < vt.nrows(); i++) {
    FAIL_IF_NOT(vt.ncols(i) == field[i].size());
    for (size_t j = 0; j < vt.ncols(i); j++)
      FAIL_IF_NOT(compare_vdf_field(vt(i, j), field[i][j]));
  }
  if (ut.numFails == 0)
    PASSMSG("test_FT passes.");
  return;
}

//------------------------------------------------------------------------------------------------//
int main(int argc, char *argv[]) {
  rtt_dsxx::ScalarUnitTest ut(argc, argv, rtt_dsxx::release);
  try {
    // >>> UNIT TESTS
    test_vector<vector<vector<int>>>(ut);
    test_vector<vector<vector<double>>>(ut);
    test_vector<vector<vector<float>>>(ut);
    test_FT<int>(ut);
    test_FT<double>(ut);
  }
  UT_EPILOG(ut);
}

//------------------------------------------------------------------------------------------------//
// end of tstViz_Traits.cc
//------------------------------------------------------------------------------------------------//
