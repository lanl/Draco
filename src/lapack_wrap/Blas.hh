//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   lapack_wrap/Blas.hh
 * \brief  Header for BLAS functions.
 * \note   Copyright (C) 2010-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#ifndef rtt_lapack_wrap_Blas_hh
#define rtt_lapack_wrap_Blas_hh

#include "Blas_Prototypes.hh"
#include "ds++/Assert.hh"
#include <algorithm>
#include <typeinfo>
#include <vector>

namespace rtt_lapack_wrap {

//------------------------------------------------------------------------------------------------//
// y <- x (COPY)
//------------------------------------------------------------------------------------------------//
/*!
 * \brief Do \f$ y\leftarrow x \f$ for type float.
 *
 * Results are written into y.
 */
inline void blas_copy(int N, const float *x, int increment_x, float *y, int increment_y) {
  Check(N >= 0);
  Check(x);
  Check(y);

  // do a single precision copy
  FC_GLOBAL(scopy, SCOPY)
  (&N, const_cast<float *>(x), &increment_x, y, &increment_y);
}

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Do \f$ y\leftarrow x \f$ for type double.
 *
 * Results are written into y.
 */
inline void blas_copy(int N, const double *x, int increment_x, double *y, int increment_y) {
  Check(N >= 0);
  Check(x);
  Check(y);

  // do a double precision axpy
  FC_GLOBAL(dcopy, DCOPY)
  (&N, const_cast<double *>(x), &increment_x, y, &increment_y);
}

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Do \f$ y\leftarrow x \f$ for vector<T> type.
 *
 * \param x vector<T> array
 * \param increment_x x stride
 * \param y vector<T> array
 * \param increment_y y stride
 *
 * The results are written into y.
 */
template <typename T>
inline void blas_copy(const std::vector<T> &x, int increment_x, std::vector<T> &y,
                      int increment_y) {
  Check(x.size() == y.size());
  Check(typeid(T) == typeid(float) || typeid(T) == typeid(double));

  blas_copy(static_cast<int>(x.size()), &x[0], increment_x, &y[0], increment_y);
}

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Do \f$ x\leftarrow\alpha x \f$ for type double.
 */
inline void blas_scal(int N, double alpha, double *x, int increment_x) {
  Check(N >= 0);
  Check(x);

  FC_GLOBAL(dscal, DSCAL)(&N, &alpha, x, &increment_x);
}

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Do \f$ x\leftarrow\alpha x \f$ for vector<T> type.
 *
 * \param x vector<T> array
 * \param alpha scalar of type T
 *
 * The results are written into x.
 */
template <typename T> inline void blas_scal(T alpha, std::vector<T> &x, int /*increment_x*/) {
  Check(typeid(T) == typeid(float) || typeid(T) == typeid(double));

  blas_scal(static_cast<int>(x.size()), alpha, &x[0], 1);
}

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Do \f$ \mbox{dot}\leftarrow x^{T}y \f$ for type double.
 */
inline double blas_dot(int N, const double *x, int increment_x, const double *y, int increment_y) {
  Check(N >= 0);
  Check(x);
  Check(y);

  // do a double precision dot (inner) product
  return FC_GLOBAL(ddot, DDOT)(&N, const_cast<double *>(x), &increment_x, const_cast<double *>(y),
                               &increment_y);
}

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Do \f$ \mbox{dot}\leftarrow x^{T}y \f$ for vector<T> type.
 *
 * \param x vector<T> array
 * \param increment_x x stride
 * \param y vector<T> array
 * \param increment_y y stride
 *
 * \return the dot product (type T)
 */
template <typename T>
inline T blas_dot(const std::vector<T> &x, int increment_x, const std::vector<T> &y,
                  int increment_y) {
  Check(x.size() == y.size());
  Check(typeid(T) == typeid(float) || typeid(T) == typeid(double));

  return blas_dot(static_cast<int>(x.size()), &x[0], increment_x, &y[0], increment_y);
}

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Do \f$ y\leftarrow\alpha x + y \f$ for type double.
 *
 * Results are written into y.
 */
inline void blas_axpy(int N, double alpha, const double *x, int increment_x, double *y,
                      int increment_y) {
  Check(N >= 0);
  Check(x);
  Check(y);

  // do a double precision axpy
  FC_GLOBAL(daxpy, DAXPY)
  (&N, &alpha, const_cast<double *>(x), &increment_x, y, &increment_y);
}

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Do \f$ y\leftarrow\alpha x + y \f$ for vector<T> type.
 *
 * \param x vector<T> array
 * \param increment_x x stride
 * \param y vector<T> array
 * \param increment_y y stride
 * \param alpha scalar of type T
 *
 * The results are written into y.
 */
template <typename T>
inline void blas_axpy(T alpha, const std::vector<T> &x, int increment_x, std::vector<T> &y,
                      int increment_y) {
  Check(x.size() == y.size());
  Check(typeid(T) == typeid(float) || typeid(T) == typeid(double));

  blas_axpy(static_cast<int>(x.size()), alpha, &x[0], increment_x, &y[0], increment_y);
}

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Do \f$ \mbox{nrm2}\leftarrow \| x\|_{2} \f$ for type double.
 */
inline double blas_nrm2(int N, const double *x, int increment_x) {
  Check(N >= 0);
  Check(x);

  // do a double precision 2-norm
  double nrm2 = FC_GLOBAL(dnrm2, DNRM2)(&N, const_cast<double *>(x), &increment_x);
  Check(nrm2 >= 0.0);
  return nrm2;
}

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Do \f$ \mbox{nrm2}\leftarrow \| x\|_{2} \f$ stl-algorithms style.
 *
 * The iterators must point to float or double value_type's.
 *
 * \param x_begin iterator pointing to the beginning of x
 * \param x_end iterator pointing to the end of x
 *
 * \return the 2-norm of x (the value_type of Forward_Iterator)
 */
template <typename Forward_Iterator>
inline typename std::iterator_traits<Forward_Iterator>::value_type
blas_nrm2(Forward_Iterator x_begin, Forward_Iterator x_end) {
  Check(typeid(typename std::iterator_traits<Forward_Iterator>::value_type) == typeid(double) ||
        typeid(typename std::iterator_traits<Forward_Iterator>::value_type) == typeid(float));

  // get the size of the arrays
  auto N = std::distance(x_begin, x_end);

  // allocate x array
  typename std::iterator_traits<Forward_Iterator>::value_type *x;
  x = new typename std::iterator_traits<Forward_Iterator>::value_type[N];

  // the dot product
  typename std::iterator_traits<Forward_Iterator>::value_type nrm2 = 0.0;

  // copy into x
  std::copy(x_begin, x_end, &x[0]);

  // do the 2-norm
  Check(N < INT_MAX);
  nrm2 = blas_nrm2(static_cast<int>(N), x, 1);
  Check(nrm2 >= 0.0);

  // clean up the memory
  delete[] x;

  // return the result
  return nrm2;
}

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Do \f$ \mbox{nrm2}\leftarrow \| x\|_{2} \f$ vector<T> type.
 *
 * \param x vector<T> array
 * \param increment_x x stride
 *
 * \return the 2-norm of x (type T)
 */
template <typename T> inline T blas_nrm2(const std::vector<T> &x, int increment_x) {
  Check(typeid(T) == typeid(float) || typeid(T) == typeid(double));

  return blas_nrm2(static_cast<int>(x.size()), &x[0], increment_x);
}

} // end namespace rtt_lapack_wrap

#endif // rtt_lapack_wrap_Blas_hh

//------------------------------------------------------------------------------------------------//
// end of lapack_wrap/Blas.hh
//------------------------------------------------------------------------------------------------//
