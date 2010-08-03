//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   traits/Viz_Traits.hh
 * \author Thomas M. Evans
 * \date   Fri Jan 21 17:10:54 2000
 * \brief  Viz_Traits header file.
 * \note   Copyright (C) 2000-2010 Los Alamos National Security, LLC.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef __traits_Viz_Traits_hh__
#define __traits_Viz_Traits_hh__

#include "ds++/Assert.hh"
#include <vector>

namespace rtt_traits
{
 
//===========================================================================//
/*!
 * \class Viz_Traits
 *
 * \brief Traits that are used by the rtt_viz package.
 *
 * These traits provide a common way to access 2D-styles arrays/fields in the
 * viz package.  Essentially, they allow rtt_dsxx::Mat2 and vector<vector<T>
 * > types to access data using (i,j) operator overloading.  There is a
 * general field templated type class; specializations exist for
 * vector<vector>.  Other specializations can be added as needed.
 *
 * The generalized class requires the Field Type (FT) template argument to
 * have the following services:
 *
 * \arg operator()(int i, int j) where the range is [0:N-1, 0:N-1]; \arg
 * *nrows() returns the number of rows (i index); \arg ncols(int row) returns
 * *the number of columns in row (j index); \arg FT::value_type defined to
 * the type returned by the field (int, double, etc).  
 */
// revision history:
// -----------------
// 0) original
// 1) 28-JAN-00 : added explicit specializations for vector<vector<int>> and
//                vector<vector<double>> because of totalview goofiness
// 
//===========================================================================//

template<class FT>
class Viz_Traits 
{
  private:
    //! Reference to the field.
    const FT &field;

  public:
    //! Constructor.
    Viz_Traits(const FT &field_in) : field(field_in) {/*...*/}

    //! Overloaded operator().
    typename FT::value_type operator()(int i, int j) const
    { return field(i, j); }

    //! Row size accessor.
    int nrows() const { return field.nrows(); }

    //! Column size accessor.
    int ncols(int row) const { return field.ncols(row); }
};

//---------------------------------------------------------------------------//
// Specialization for std::vector<std::vector>

template<class T>
class Viz_Traits< std::vector<std::vector<T> > >
{
  private:
    // Reference to vector<vector> field.
    const std::vector<std::vector<T> > &field;
    
  public:
    // Constructor.
    Viz_Traits(const std::vector<std::vector<T> > &fin) : field(fin)
    {
	// Nothing to do here
    } 

    // Overloaded operator().
    T operator()(int i, int j) const
    {
	Require(i >= 0 && i < field.size());
	Require(j >= 0 && j < field[i].size());
	return field[i][j];
    }

    // Row size accessor.
    int nrows() const { return field.size(); }

    // Column size accessor.
    int ncols(int row) const
    {
	Require (row >= 0 && row < field.size());
	return field[row].size();
    }
};

//---------------------------------------------------------------------------//
// explicit specialization for vector<vector<int>> --> Need this because of
// totalviews inability to handle templated specializations

template<>
class Viz_Traits< std::vector<std::vector<int> > >
{
  private:
    // Reference to vector<vector> field.
    const std::vector<std::vector<int> > &field;
    
  public:
    // Constructor.
    Viz_Traits(const std::vector<std::vector<int> > &fin) : field(fin)
    {
	// Nothing to do here
    } 

    // Overloaded operator().
    int operator()(int i, int j) const
    {
	Require(i >= 0 && i < static_cast<int>(field.size()));
	Require(j >= 0 && j < static_cast<int>(field[i].size()));
	return field[i][j];
    }

    // Row size accessor.
    int nrows() const { return field.size(); }

    // Column size accessor.
    int ncols(int row) const
    {
	Require (row >= 0 && row < static_cast<int>(field.size()));
	return field[row].size();
    }
};

//---------------------------------------------------------------------------//
// explicit specialization for vector<vector<double>> --> Need this because of
// totalviews inability to handle templated specializations

template<>
class Viz_Traits< std::vector<std::vector<double> > >
{
  private:
    // Reference to vector<vector> field.
    const std::vector<std::vector<double> > &field;
    
  public:
    // Constructor.
    Viz_Traits(const std::vector<std::vector<double> > &fin) : field(fin)
    {
	// Nothing to do here
    } 

    // Overloaded operator().
    double operator()(int i, int j) const
    {
	Require(i >= 0 && i < static_cast<int>(field.size()));
	Require(j >= 0 && j < static_cast<int>(field[i].size()));
	return field[i][j];
    }

    // Row size accessor.
    int nrows() const { return field.size(); }

    // Column size accessor.
    int ncols(int row) const
    {
	Require (row >= 0 && row < static_cast<int>(field.size()));
	return field[row].size();
    }
};

} // end namespace rtt_traits

#endif                          // __traits_Viz_Traits_hh__

//---------------------------------------------------------------------------//
//                              end of traits/Viz_Traits.hh
//---------------------------------------------------------------------------//
