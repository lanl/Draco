//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   c4/test/tstSLURM_Task_Info.cc
 * \author Tim Kelley
 * \date   Fri Jun 7 08:06:53 2019
 * \note   Copyright (C) 2020-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "c4/SLURM_Task_Info.hh"
#include "ds++/Release.hh"
#include "ds++/ScalarUnitTest.hh"
#include "ds++/SystemCall.hh"
#include <array>
#include <cstdlib>
#include <functional> // std::function
#include <map>
#include <thread>

using rtt_dsxx::UnitTest;

using env_store_value = std::pair<bool, std::string>;
using env_store_t = std::map<std::string, env_store_value>;

//------------------------------------------------------------------------------------------------//
/* Helper function: Record SLURM keys and values, if any, then remove them from environment. Return
 * recorded values so they can be restored later. */
env_store_t clean_env() {
  // for each key, is it defined? If so, record the value, then unset it. If not, note that.
  std::array<std::string, 3> slurm_keys = {"SLURM_CPUS_PER_TASK", "SLURM_NTASKS",
                                           "SLURM_JOB_NUM_NODES"};
  env_store_t store{};
  for (auto &k : slurm_keys) {
    // We're assuming none of those flags were set to something like 2 billion
    constexpr int a_large_int = 0x0FFFFFFC;
    bool k_defined{false};
    int ival;
    std::tie(k_defined, ival) = rtt_dsxx::get_env_val<int>(k, a_large_int);
    if (k_defined) {
      std::stringstream storestr;
      storestr << ival;
      store.insert({k, {true, storestr.str()}});
      int unset_ok = draco_unsetenv(k.c_str());
      if (0 != unset_ok) {
        std::cout << __FUNCTION__ // NOLINT
                  << ":" << __LINE__ << "Failed to unset environment variable " << k
                  << "! errno = " << errno << std::endl;
      }
    } else {
      store.insert({k, {false, ""}});
    }
  } // for k in slurm keys
  return store;
} // clean_env

//------------------------------------------------------------------------------------------------//
/* Helper function: Restore the SLURM keys that were previously defined. */
void restore_env(env_store_t const &store) {
  /* For each key, if it was defined, restore that definition. If it was not defined, destroy any
   * subsequent definition */
  for (auto const &str_it : store) {
    std::string const &key{str_it.first};
    env_store_value const &val{str_it.second};
    bool const &was_defined{val.first};
    if (was_defined) {
      std::string const &val_str{val.second};
      int set_ok = draco_setenv(key.c_str(), val_str.c_str());
      if (0 != set_ok) {
        std::cout << __FUNCTION__ // NOLINT
                  << ":" << __LINE__ << "Failed to set environment variable! errno = " << errno
                  << std::endl;
      }
    } else {
      int unset_ok = draco_unsetenv(key.c_str());
      if (0 != unset_ok) {
        std::cout << __FUNCTION__ // NOLINT
                  << ":" << __LINE__ << "Failed to unset environment variable! errno = " << errno
                  << std::endl;
      }
    }
  } // for things in store
  return;
} // restore_env

//------------------------------------------------------------------------------------------------//
/* Test with a "clean" environment--that is, no slurm keys. */
void test_instantiate_SLURM_Info(UnitTest &ut) {
  /* Test instantiating SLURM_Task_Info in a 'clean' environment */
  auto env_tmp = clean_env();
  rtt_c4::SLURM_Task_Info ti;
  FAIL_IF(ti.is_cpus_per_task_set());
  FAIL_IF(ti.is_ntasks_set());
  FAIL_IF(ti.is_job_num_nodes_set());
  FAIL_IF_NOT(ti.get_cpus_per_task() == 0xFFFFFFF);
  FAIL_IF_NOT(ti.get_ntasks() == 0xFFFFFFE);
  FAIL_IF_NOT(ti.get_job_num_nodes() == 0xFFFFFFD);
  restore_env(env_tmp);
  return;
}

//------------------------------------------------------------------------------------------------//
/* Test with a "live" environment--that is, slurm keys are defined. */
void test_SLURM_Info(UnitTest &ut) {
  /* Test instantiating SLURM_Task_Info in a 'clean' environment */
  auto orig_env = clean_env();
  int const iset_cpus_per_task{21}, iset_ntasks{341}, iset_job_num_nodes{1001};
  const char *set_cpus_per_task{"21"}, *set_ntasks{"341"}, *set_job_num_nodes{"1001"};
  env_store_t test_env = {{"SLURM_CPUS_PER_TASK", {true, set_cpus_per_task}},
                          {"SLURM_NTASKS", {true, set_ntasks}},
                          {"SLURM_JOB_NUM_NODES", {true, set_job_num_nodes}}};
  restore_env(test_env);
  rtt_c4::SLURM_Task_Info ti;
  FAIL_IF_NOT(ti.is_cpus_per_task_set());
  FAIL_IF_NOT(ti.is_ntasks_set());
  FAIL_IF_NOT(ti.is_job_num_nodes_set());
  FAIL_IF_NOT(ti.get_cpus_per_task() == iset_cpus_per_task);
  FAIL_IF_NOT(ti.get_ntasks() == iset_ntasks);
  FAIL_IF_NOT(ti.get_job_num_nodes() == iset_job_num_nodes);
  restore_env(orig_env);
  return;
}

//------------------------------------------------------------------------------------------------//
/* Test with a partial "live" environment--that is, slurm keys are defined. */
void test_SLURM_Info_partial(UnitTest &ut) {
  /* Test instantiating SLURM_Task_Info in a 'clean' environment */
  auto orig_env = clean_env();
  int const iset_cpus_per_task{21}, iset_job_num_nodes{1001};
  const char *set_cpus_per_task{"21"}, *set_job_num_nodes{"1001"};
  env_store_t test_env = {{"SLURM_CPUS_PER_TASK", {true, set_cpus_per_task}},
                          {"SLURM_JOB_NUM_NODES", {true, set_job_num_nodes}}};
  restore_env(test_env);
  rtt_c4::SLURM_Task_Info ti;
  FAIL_IF_NOT(ti.is_cpus_per_task_set());
  FAIL_IF(ti.is_ntasks_set());
  FAIL_IF_NOT(ti.is_job_num_nodes_set());
  FAIL_IF_NOT(ti.get_cpus_per_task() == iset_cpus_per_task);
  FAIL_IF_NOT(ti.get_ntasks() == 0xFFFFFFE);
  FAIL_IF_NOT(ti.get_job_num_nodes() == iset_job_num_nodes);
  restore_env(orig_env);
  return;
}

using t_func = std::function<void(UnitTest &)>;

//------------------------------------------------------------------------------------------------//
void run_a_test(UnitTest &u, t_func const &f, std::string const &msg) {
  f(u);
  if (u.numFails == 0) {
    u.passes(msg);
  }
  return;
}

//------------------------------------------------------------------------------------------------//
int main(int argc, char *argv[]) {
  rtt_dsxx::ScalarUnitTest ut(argc, argv, rtt_dsxx::release);
  try {
    run_a_test(ut, test_instantiate_SLURM_Info, "SLURM_Info (clean env) ok.");
    run_a_test(ut, test_SLURM_Info, "SLURM_Info (SLURM vars set) ok.");
    run_a_test(ut, test_SLURM_Info_partial, "SLURM_Info (partial SLURM vars set) ok.");
  } // try--catches in the epilog:
  UT_EPILOG(ut);
}

//------------------------------------------------------------------------------------------------//
// end of c4/test/tstSLURM_Task_Info.cc
//------------------------------------------------------------------------------------------------//
