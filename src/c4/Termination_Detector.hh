//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   c4/Termination_Detector.hh
 * \author Kent Budge
 * \brief  Definition of class Termination_Detector
 * \note   Copyright (C) 2010-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#ifndef c4_Termination_Detector_hh
#define c4_Termination_Detector_hh

#include "ds++/config.h"

namespace rtt_c4 {

//================================================================================================//
/*!
 * \class Termination_Detector
 * \brief Detect termination of an indeterminate parallel algorithm.
 *
 * This class is our implementation in C++ of Will McLendon's (Sandia) termination package.  It is
 * used to detect termination of an indeterminate algorithm.
 *
 * An indeterminate algorithm is one for which it is not possible to determine in advance how much
 * work must be done on each processor.  One example of an indeterminate algorithm is traversal of a
 * directed graph that may have cycles. A processor performing such a sweep cannot know in advance
 * whether all of its nodes will be visited by the traversal, which makes the traversal
 * indeterminate.
 *
 * The code is used as follows: Every time a processor is able to make progress on the computation
 * (as measured in arbitrary units of work) the Update_Work_Count function should be called to
 * indicate the progress made. Whenever data is communicated to another processor, the
 * Update_Send_Count function should be called to indicate the messages sent out. And whenever data
 * is received, the Update_Receive_Count function should be called to indicate the messages
 * received.
 *
 * Whenever a processor has no more work to do, it should be sure the counts are all updated, then
 * call Process() to see if the algorithm has terminated. Termination occurs when no processor has
 * work to do (as implied by the call to Process()), no processor has done any work since the last
 * check, and all sent messages have been received.
 *
 * See test/tstTermination_Detector for an example of how this works.
 */
//================================================================================================//

class Termination_Detector {
public:
  // NESTED CLASSES AND TYPEDEFS

  // CREATORS

  //! Constructor
  explicit Termination_Detector(int tag);

  //! Destructor.
  ~Termination_Detector(); // NOLINT(performance-trivially-destructible) bc optional debug behavior.

  //! Copy/Move assignment/construction operators not implemented.
  Termination_Detector &operator=(const Termination_Detector &rhs) = delete;
  Termination_Detector &operator=(const Termination_Detector &&rhs) = delete;
  Termination_Detector(const Termination_Detector &rhs) = delete;
  Termination_Detector(const Termination_Detector &&rhs) = delete;

  // MANIPULATORS

  //! Indicate that an indetermine algorithm is starting.
  void init();

  //! Indicate that a certain number of units of work have been performed.
  void update_work_count(unsigned units_of_work) { work_count_ += units_of_work; }

  //! Indicate that a certain number of messages have been received.
  void update_receive_count(unsigned messages_received) { receive_count_ += messages_received; }

  //! Indicate that a certain number of messages have been sent.
  void update_send_count(unsigned messages_sent) { send_count_ += messages_sent; }

  //! See if the algorithm has terminated.
  bool is_terminated();

  // ACCESSORS

private:
  // NESTED CLASSES AND TYPEDEFS

  //! What is the state of this processor?
  enum State { DOWN, UP, TERMINATED };

  //! What action is the parent of this processor requesting?
  enum Parent_Action_Request { SEND_DATA, TERMINATE };

  //! What sort of processor is this?
  enum Processor_Type { ROOT, LEAF, INTERNAL };

  // DATA

  int tag_;
  unsigned number_of_processors_;
  unsigned pid_;
  unsigned parent_pid_, son_pid_, daughter_pid_;
  Processor_Type ptype_;
  State state_;
  unsigned send_count_, receive_count_, work_count_;
  unsigned subtree_send_count_, subtree_receive_count_, subtree_work_count_;
  unsigned old_global_work_count_;
};

} // end namespace rtt_c4

#endif // c4_Termination_Detector_hh

//------------------------------------------------------------------------------------------------//
// end of c4/Termination_Detector.hh
//------------------------------------------------------------------------------------------------//
