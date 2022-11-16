//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   ds++/SortPermutation.hh
 * \author Randy M. Roberts
 * \date   Mon Feb 14 14:18:27 2000
 * \brief  SortPermutation class definition.
 * \note   Copyright (C) 2010-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#ifndef rtt_dsxx_SortPermutation_hh
#define rtt_dsxx_SortPermutation_hh

#include "Assert.hh"
#include "isSorted.hh"
#include <algorithm>
#include <functional>
#include <iterator>
#include <vector>

namespace rtt_dsxx {

//================================================================================================//
/*!
 * \class SortPermutation
 *
 * This class determines a permutation used to sort a sequence.
 *
 * This is necessary if multiple sequences need to be re-arranged, via the ordering necessary to
 * sort a sequece of "keys".
 *
 * Create an object of class SortPermutation using either of the following constructors:
 * \code
 *      template<typename ForwardIterator>
 *      SortPermutation(ForwardIterator first, ForwardIterator last)
 *
 *      template<typename ForwardIterator, class StrictWeakOrdering>
 *      SortPermutation(ForwardIterator first, ForwardIterator last,
 *                      StrictWeakOrdering comp)
 * \endcode
 * The object creates a permutation that results in an ordering of elements in [first,last) into
 * nondescending order.  This ordering is not necessarily stable.
 *
 * The SortPermutation object can be used to access the "index table" for the sequence via
 * operator[](int), begin(), and end().
 *
 * An index table is a table of integer pointers telling which number sequence element comes first
 * in ascending order, which second, and so on.  See "Numerical Recipes" for a full discussion of an
 * index table.
 *
 * The inverse operations on the SortPermutation object can be used to access the "rank table" for
 * the sequence via inv(int), inv_begin(), and inv_end().
 *
 * A rank table is a table telling what the numerical rank of the first sequence element, the second
 * sequence element, and so on.  See "Numerical Recipes" for a full discussion of a rank table.
 */
/*!
 * \example tstSortPermutation.cc
 *
 * Test of rtt_dsxx::SortPermutation and isSorted.hh functions.
 */
//================================================================================================//

class SortPermutation {

  // NESTED CLASSES AND TYPEDEFS

public:
  typedef int value_type; // NOLINT

private:
  typedef std::vector<value_type> InternalRep; // NOLINT

public:
  // A SortPermutation can not be modified; therefor, always use a
  // const_iterator.

  typedef InternalRep::const_iterator iterator;       // NOLINT
  typedef InternalRep::const_iterator const_iterator; // NOLINT

private:
  // Forward Declarations

  template <typename COMP> class CompareProxy;
  template <typename IT> class Proxy;

  template <typename COMP> friend class CompareProxy;
  template <typename IT> friend class Proxy;

  template <typename IT> class Proxy {
    friend class CompareProxy<Proxy>;

    typedef typename std::iterator_traits<IT>::value_type value_type; // NOLINT

    SortPermutation::value_type pos;
    const std::vector<IT> &iters;

  public:
    Proxy(SortPermutation::value_type pos_, const std::vector<IT> &iters_)
        : pos(pos_), iters(iters_) {}
    ~Proxy() = default;
    Proxy(Proxy const &rhs) noexcept : pos(rhs.pos), iters(rhs.iters) {}
    Proxy(Proxy &&rhs) noexcept : pos(rhs.pos), iters(rhs.iters) {}
    Proxy &operator=(Proxy const &rhs) {
      pos = rhs.pos;
      return *this;
    }
    Proxy &operator=(Proxy &&rhs) noexcept {
      pos = rhs.pos;
      return *this;
    }

    const value_type &value() const { return *iters[pos]; }

    operator SortPermutation::value_type() { return pos; } // NOLINT [hicpp-explicit-conversions]
  };

  template <typename COMP> class CompareProxy {
  public:
    const COMP &comp;
    explicit CompareProxy(const COMP &comp_) : comp(comp_) {}
    template <typename IT> CompareProxy &operator=(CompareProxy const &comp_);
    template <typename IT> bool operator()(const Proxy<IT> &p1, const Proxy<IT> &p2) const {
      return comp(p1.value(), p2.value());
    }
  };

  // DATA

private:
  InternalRep indexTable_m;
  InternalRep rankTable_m;

public:
  // CREATORS

  template <typename IT, class COMP>
  SortPermutation(IT first, IT last, const COMP &comp)
      : indexTable_m(std::distance(first, last)), rankTable_m(indexTable_m.size()) {
    createPermutation(first, last, comp);
  }

  template <typename IT>
  SortPermutation(IT first, IT last)
      : indexTable_m(std::distance(first, last)), rankTable_m(indexTable_m.size()) {
    typedef typename std::iterator_traits<IT>::value_type vtype; // NOLINT
    createPermutation(first, last, std::less<vtype>());
  }

  // ACCESSORS

  /*!
   * \brief Returns the i'th entry into the index table.
   * \param i The i'th entry into the sorted order.  The condition, *(first + sortPerm[i+1]) <
   *        *(first + sortPerm[i]), is guaranteed to be false.
   *
   * For example,
   * \code
   *     first = unsorted.begin();
   *     last  = unsorted.end()
   *     SortPermutation sortPerm(first, last);
   *     for (int i=0; i<unsorted.size(); i++)
   *        sorted[i] = unsorted[sortPerm[i]];
   * \endcode
   * results in sorted containing the sorted elements of [first, last).
   */

  value_type operator[](unsigned i) const { return indexTable_m[i]; }

  //! Returns the begin const_iterator into the index table.
  const_iterator begin() const { return indexTable_m.begin(); }

  //! Returns the end const_iterator into the index table.
  const_iterator end() const { return indexTable_m.end(); }

  /*!
   * \brief Returns the i'th entry into the rank table.
   * \param i The i'th entry from the sorted order.
   *
   * For example,
   * \code
   *     first = unsorted.begin();
   *     last  = unsorted.end()
   *     SortPermutation sortPerm(first, last);
   *     for (int i=0; i<unsorted.size(); i++)
   *        sorted[sortPerm.inv(i)] = unsorted[i];
   * \endcode
   * results in sorted containing the sorted elements of [first, last).
   */
  value_type inv(int i) const { return rankTable_m[i]; }

  //! Returns the begin const_iterator into the rank table.
  const_iterator inv_begin() const { return rankTable_m.begin(); }

  //! Returns the end const_iterator into the rank table.
  const_iterator inv_end() const { return rankTable_m.end(); }

  //! Returns the size of the index and rank tables.
  int size() const { return static_cast<int>(indexTable_m.size()); }

private:
  // IMPLEMENTATION

  template <typename IT, class COMP> void createPermutation(IT first, IT last, const COMP &comp) {
    std::vector<IT> iters;
    iters.reserve(size());
    IT it = first;
    while (it != last) {
      iters.push_back(it);
      ++it;
    }
    doCreatePermutation(first, last, comp, iters);
  }

#ifdef ENSURE_ON
  template <typename IT, class COMP> bool isPermutationSorted(IT first, IT last, const COMP &comp) {
    typedef typename std::iterator_traits<IT>::value_type vtype; // NOLINT
    std::vector<vtype> vv(first, last);

    for (int i = 0; first != last && i < size(); ++i, ++first) {
      vv[inv(i)] = *first;
    }

    return isSorted(vv.begin(), vv.end(), comp);
  }
#endif

  template <typename IT, class COMP>
  void doCreatePermutation(IT Remember(first), IT Remember(last), const COMP &comp,
                           const std::vector<IT> &iters) {
    std::vector<Proxy<IT>> proxies;
    proxies.reserve(size());

    for (SortPermutation::value_type i = 0; i < size(); ++i)
      proxies.push_back(Proxy<IT>(i, iters));

    std::sort(proxies.begin(), proxies.end(), CompareProxy<COMP>(comp));

    for (SortPermutation::value_type i = 0; i < size(); ++i) {
      indexTable_m[i] = proxies[i];
      rankTable_m[indexTable_m[i]] = i;
    }

    Ensure(isPermutationSorted(first, last, comp));
  }
};

} // end namespace rtt_dsxx

#endif // rtt_dsxx_SortPermutation_hh

//------------------------------------------------------------------------------------------------//
// end of ds++/SortPermutation.hh
//------------------------------------------------------------------------------------------------//
