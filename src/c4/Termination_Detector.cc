//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   c4/Termination_Detector.cc
 * \author Kent Budge
 * \date   Thu Jan 12 10:27:45 2006
 * \note   Copyright (C) 2010-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "Termination_Detector.hh"
#include "C4_Functions.hh"
#include <array>

#undef PRINTF_DEBUG

namespace rtt_c4 {

//------------------------------------------------------------------------------------------------//
/*
 * \brief Constructor
 * \param[in] tag Tag to use for termination algorithm messages. This must be distinct from any
 *                other message tags that may simultaneously be in use.
 */
Termination_Detector::Termination_Detector(int const tag)
    : tag_(tag), number_of_processors_(nodes()), pid_(node()), parent_pid_((pid_ - 1) / 2),
      son_pid_(2 * pid_ + 1), daughter_pid_(son_pid_ + 1), ptype_(), state_(UP), send_count_(0),
      receive_count_(0), work_count_(0), subtree_send_count_(0), subtree_receive_count_(0),
      subtree_work_count_(0), old_global_work_count_(0) {
  if (pid_ == 0) {
    ptype_ = ROOT;
  } else {
    if (son_pid_ >= number_of_processors_ && daughter_pid_ >= number_of_processors_) {
      ptype_ = LEAF;
    } else {
      ptype_ = INTERNAL;
    }
  }
#ifdef PRINTF_DEBUG
  cout << pid_ << " constructed with tag " << tag_ << endl;
#endif
}

//------------------------------------------------------------------------------------------------//
void Termination_Detector::init() {
  state_ = UP;
  send_count_ = 0;
  receive_count_ = 0;
  work_count_ = 0;
  old_global_work_count_ = 0;
}

//------------------------------------------------------------------------------------------------//
#ifdef PRINTF_DEBUG
Termination_Detector::~Termination_Detector() {
  cout << pid_ << '/' << tag_ << ": destroyed" << endl;
}
#else
Termination_Detector::~Termination_Detector() = default;
#endif

//------------------------------------------------------------------------------------------------//
//! Return \c true if we have terminated; \c false otherwise.
bool Termination_Detector::is_terminated() {

  // We keep track of three integral values in this buffer container.
  // - send_count
  // - receive_count
  // - work_count
  size_t constexpr num_counts = 3;
  std::array<unsigned, num_counts> buffer{0, 0, 0};

  if (ptype_ == ROOT)
  // root processor
  {
    if (son_pid_ < number_of_processors_) {
#ifdef PRINTF_DEBUG
      cout << "0/" << tag_ << ": expecting from " << son_pid_ << endl;
#endif
      receive(buffer.data(), num_counts, son_pid_, tag_);
#ifdef PRINTF_DEBUG
      cout << "0/" << tag_ << ": received " << son_pid_ << ": " << buffer[0] << ' ' << buffer[1]
           << ' ' << buffer[2] << endl;
#endif

      subtree_send_count_ = buffer[0];
      subtree_receive_count_ = buffer[1];
      subtree_work_count_ = buffer[2];
    }
    if (daughter_pid_ < number_of_processors_) {
#ifdef PRINTF_DEBUG
      cout << "0/" << tag_ << ": expecting from " << daughter_pid_ << endl;
#endif
      receive(buffer.data(), num_counts, daughter_pid_, tag_);
#ifdef PRINTF_DEBUG
      cout << "0/" << tag_ << ": received " << daughter_pid_ << ": " << buffer[0] << ' '
           << buffer[1] << ' ' << buffer[2] << endl;
#endif

      // If the daughter exists, the son also exists, and the subtree counts are already
      // initialized.
      subtree_send_count_ += buffer[0];
      subtree_receive_count_ += buffer[1];
      subtree_work_count_ += buffer[2];
    }

    unsigned global_send_count = subtree_send_count_ + send_count_;

#ifdef PRINTF_DEBUG
    cout << "0/" << tag_ << ": global send count: " << global_send_count << endl;
#endif

    unsigned global_receive_count = subtree_receive_count_ + receive_count_;

#ifdef PRINTF_DEBUG
    cout << "0/" << tag_ << ": global receive count: " << global_receive_count << endl;
#endif

    unsigned global_work_count = subtree_work_count_ + work_count_;

#ifdef PRINTF_DEBUG
    cout << "0/" << tag_ << ": global work count: " << global_work_count << endl;
#endif

    buffer[0] =
        (global_work_count == old_global_work_count_ && global_send_count == global_receive_count)
            ? TERMINATE
            : SEND_DATA;

    old_global_work_count_ = global_work_count;

    if (son_pid_ < number_of_processors_) {
      send(buffer.data(), 1, son_pid_, tag_);

#ifdef PRINTF_DEBUG
      cout << "0/" << tag_ << ": sent " << son_pid_ << ": " << buffer[0] << endl;
#endif
    }
    if (daughter_pid_ < number_of_processors_) {
      send(buffer.data(), 1, daughter_pid_, tag_);

#ifdef PRINTF_DEBUG
      cout << "0/" << tag_ << ": sent " << daughter_pid_ << ": " << buffer[0] << endl;
#endif
    }
    if (buffer[0] == TERMINATE) {
#ifdef PRINTF_DEBUG
      cout << pid_ << '/' << tag_ << ": destroyed" << endl;
#endif
      return true;
    }
  } else if (ptype_ == LEAF) {
    if (state_ == UP) {
      buffer[0] = send_count_;
      buffer[1] = receive_count_;
      buffer[2] = work_count_;

      send(buffer.data(), num_counts, parent_pid_, tag_);

      state_ = DOWN;

#ifdef PRINTF_DEBUG
      cout << pid_ << '/' << tag_ << ": sent " << parent_pid_ << ": " << buffer[0] << ' '
           << buffer[1] << ' ' << buffer[2] << ", state now DOWN" << endl;
#endif
    } else {
#ifdef PRINTF_DEBUG
      cout << pid_ << '/' << tag_ << ": expecting from " << parent_pid_ << endl;
#endif
      receive(buffer.data(), 1, parent_pid_, tag_);
#ifdef PRINTF_DEBUG
      cout << pid_ << '/' << tag_ << ": received " << parent_pid_ << ": " << buffer[0] << endl;
#endif

      if (buffer[0] == TERMINATE) {
#ifdef PRINTF_DEBUG
        cout << pid_ << '/' << tag_ << ": destroyed" << endl;
#endif
        return true;
      } else {
        buffer[0] = send_count_;
        buffer[1] = receive_count_;
        buffer[2] = work_count_;

        send(buffer.data(), num_counts, parent_pid_, tag_);

#ifdef PRINTF_DEBUG
        cout << pid_ << '/' << tag_ << ": sent " << parent_pid_ << ": " << buffer[0] << ' '
             << buffer[1] << ' ' << buffer[2] << ", state still DOWN" << endl;
#endif
      }
    }
  } else
  // internal processor
  {
    if (state_ == UP) {
      if (son_pid_ < number_of_processors_) {
#ifdef PRINTF_DEBUG
        cout << pid_ << '/' << tag_ << ": expecting from " << son_pid_ << endl;
#endif
        receive(buffer.data(), num_counts, son_pid_, tag_);
#ifdef PRINTF_DEBUG
        cout << pid_ << '/' << tag_ << ": received " << son_pid_ << ": " << buffer[0] << ' '
             << buffer[1] << ' ' << buffer[2] << endl;
#endif

        subtree_send_count_ = buffer[0];
        subtree_receive_count_ = buffer[1];
        subtree_work_count_ = buffer[2];
      }
      if (daughter_pid_ < number_of_processors_) {
#ifdef PRINTF_DEBUG
        cout << pid_ << '/' << tag_ << ": expecting from " << daughter_pid_ << endl;
#endif
        receive(buffer.data(), num_counts, daughter_pid_, tag_);
#ifdef PRINTF_DEBUG
        cout << pid_ << '/' << tag_ << ": received " << daughter_pid_ << ": " << buffer[0] << ' '
             << buffer[1] << ' ' << buffer[2] << endl;
#endif

        // If the daughter exists, the son also exists, and the subtree counts
        // are already initialized.
        subtree_send_count_ += buffer[0];
        subtree_receive_count_ += buffer[1];
        subtree_work_count_ += buffer[2];
      }

      buffer[0] = send_count_ + subtree_send_count_;
      buffer[1] = receive_count_ + subtree_receive_count_;
      buffer[2] = work_count_ + subtree_work_count_;

      send(buffer.data(), num_counts, parent_pid_, tag_);

#ifdef PRINTF_DEBUG
      cout << pid_ << '/' << tag_ << ": sent " << parent_pid_ << ": " << buffer[0] << ' '
           << buffer[1] << ' ' << buffer[2] << ", state now DOWN" << endl;
#endif

      state_ = DOWN;
    } else
    // state==DOWN
    {
#ifdef PRINTF_DEBUG
      cout << pid_ << '/' << tag_ << ": expecting from " << parent_pid_ << endl;
#endif
      receive(buffer.data(), 1, parent_pid_, tag_);
#ifdef PRINTF_DEBUG
      cout << pid_ << '/' << tag_ << ": received " << parent_pid_ << ": " << buffer[0] << endl;
#endif

      if (son_pid_ < number_of_processors_) {
        send(buffer.data(), 1, son_pid_, tag_);

#ifdef PRINTF_DEBUG
        cout << pid_ << '/' << tag_ << ": sent " << son_pid_ << ": " << buffer[0] << endl;
#endif
      }
      if (daughter_pid_ < number_of_processors_) {
        send(buffer.data(), 1, daughter_pid_, tag_);

#ifdef PRINTF_DEBUG
        cout << pid_ << '/' << tag_ << ": sent " << daughter_pid_ << ": " << buffer[0] << endl;
#endif
      }
      state_ = UP;

#ifdef PRINTF_DEBUG
      cout << pid_ << '/' << tag_ << ": state now UP" << endl;
#endif

      if (buffer[0] == TERMINATE) {
#ifdef PRINTF_DEBUG
        cout << pid_ << '/' << tag_ << ": destroyed" << endl;
#endif
        return true;
      }
    }
  }
  return false;
}

} // end namespace rtt_c4

//------------------------------------------------------------------------------------------------//
// end of Termination_Detector.cc
//------------------------------------------------------------------------------------------------//
