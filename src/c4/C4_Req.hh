//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   c4/C4_Req.hh
 * \author Thomas M. Evans, Geoffrey Furnish
 * \date   Thu Jun  2 09:54:02 2005
 * \brief  C4_Req class definition.
 * \note   Copyright (C) 2010-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#ifndef c4_C4_Req_hh
#define c4_C4_Req_hh

// C4 package configure
#include "C4_Status.hh"
#include "C4_Traits.hh"
#include "c4/config.h"

namespace rtt_c4 {
//================================================================================================//
/*!
 * \class C4_ReqRefRep
 * \brief Handle for non-blocking message requests.
 *
 * This class provides an encapsulation for the message requests (MPI) which are produced by non
 * blocking calls.  This class automatically waits for the message to complete when the containing
 * object goes out of scope, thus plugging one of the easiest types of programming errors with non
 * blocking messaging.  Reference counting is used so that these may be passed by value without
 * accidentally triggering a program stall.
 *
 * This is a "work" class. The user interface for requests is provided by rtt_c4::C4_Req.
 */
//================================================================================================//

class C4_ReqRefRep {
  friend class C4_Req;

  // number of ref counts
  int n{0};

  // if true, we hold a request
  bool assigned{false};

#ifdef C4_MPI
  MPI_Request r;
#endif

private:
  // Private default constructor and destructor for access from C4_Req only.
  C4_ReqRefRep();

  /*! \brief Destructor.
   *
   * It is important that all existing requests are cleared before the destructor is called.  We
   * used to have a wait() in here; however, this causes exception safety problems.  In any case, it
   * is probably a bad idea to clean up communication by going out of scope. */
  ~C4_ReqRefRep() = default;

public:
  // Disallowed assignment, move assignment, copy ctor and move ctor methods.
  C4_ReqRefRep(const C4_ReqRefRep &rep) = delete;
  C4_ReqRefRep(C4_ReqRefRep &&rep) noexcept = delete;
  C4_ReqRefRep &operator=(const C4_ReqRefRep &rep) = delete;
  C4_ReqRefRep &operator=(C4_ReqRefRep &&rep) noexcept = delete;

  void wait(C4_Status *status = nullptr);
  bool complete(C4_Status *status = nullptr);
  void free();

  bool inuse() const {
#ifdef C4_MPI
    if (assigned) {
      Check(r != MPI_REQUEST_NULL);
    }
#endif
    return assigned;
  }

private:
  void set() { assigned = true; }
  void clear() { assigned = false; }
};

//================================================================================================//
/*!
 * \class C4_Req
 * \brief Non-blocking communication request class.
 *
 * This class provides an encapsulation for the message requests (MPI) which are produced by non
 * blocking calls.  This class automatically waits for the message to complete when the containing
 * object goes out of scope, thus plugging one of the easiest types of programming errors with non
 * blocking messaging.  Reference counting is used so that these may be passed by value without
 * accidentally triggering a program stall.
 *
 * This class provides an interface for non-blocking request handles that should be used by users.
 */
//================================================================================================//

class C4_Req {
  //! Request handle.
  C4_ReqRefRep *p;

public:
  C4_Req();
  C4_Req(const C4_Req &req);
  C4_Req(C4_Req &&req) noexcept;
  ~C4_Req();
  C4_Req &operator=(const C4_Req &req);
  C4_Req &operator=(C4_Req &&req) noexcept;

  //! \brief Equivalence operator
  bool operator==(const C4_Req &right) { return (p == right.p); }
  bool operator!=(const C4_Req &right) { return (p != right.p); }

  void wait(C4_Status *status = nullptr) { p->wait(status); }
  bool complete(C4_Status *status = nullptr) { return p->complete(status); }
  void free() { p->free(); }
  bool inuse() const {
#ifndef __clang_analyzer__
    Insist(p != nullptr, "attempting to use freed memory.");
    return p->inuse();
#else
    return false;
#endif
  }

private:
  void set() { p->set(); }

  // Private access to the C4_ReqRefRep internals.

#ifdef C4_MPI
  MPI_Request &r() { return p->r; }
#endif

  void free_();

  /* FRIENDSHIP
   *
   * Specific friend C4 functions that may need to manipulate the C4_ReqRefRep internals. A friend
   * function of a class is defined outside that class' scope but it has the right to access all
   * private and protected members of the class. Even though the prototypes for friend functions
   * appear in the class definition, friends are not member functions.  \ref
   * https://www.tutorialspoint.com/cplusplus/cpp_friend_functions.htm
   */
  template <typename T> friend C4_Req send_async(const T *buf, int nels, int dest, int tag);
  template <typename T>
  friend void send_async(C4_Req &r, const T *buf, int nels, int dest, int tag);
  template <typename T> friend C4_Req receive_async(T *buf, int nels, int source, int tag);
  template <typename T> friend void receive_async(C4_Req &r, T *buf, int nels, int source, int tag);
#ifdef C4_MPI
  template <typename T>
  friend void send_is_custom(C4_Req &request, T const *buffer, int size, int destination, int tag);
  template <typename T>
  friend void receive_async_custom(C4_Req &request, T *buffer, int size, int source, int tag);
  template <typename T> friend void send_is(C4_Req &r, const T *buf, int nels, int dest, int tag);
  friend void wait_all(unsigned count, C4_Req *requests);
  friend std::vector<int> wait_all_with_source(unsigned count, C4_Req *requests);
  friend unsigned wait_any(unsigned count, C4_Req *requests);
  template <typename T> friend void global_isum(T &send_buffer, T &recv_buffer, C4_Req &request);

#endif
};

} // end namespace rtt_c4

#endif // c4_C4_Req_hh

//------------------------------------------------------------------------------------------------//
// end of c4/C4_Req.hh
//------------------------------------------------------------------------------------------------//
