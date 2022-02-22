//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   ds++/test/tstXGetopt.cc
 * \author Kelly Thompson <kgt@lanl.gov>
 * \date   Wed Nov 10 09:35:09 2010
 * \brief  Test functions defined in ds++/XGetopt.cc
 * \note   Copyright (C) 2015-2022 Triad National Security, LLC., All rights reserved.
 */
//------------------------------------------------------------------------------------------------//

#include "ds++/DracoStrings.hh"
#include "ds++/Release.hh"
#include "ds++/ScalarUnitTest.hh"
#include "ds++/XGetopt.hh"
#include <array>
#include <sstream>

//------------------------------------------------------------------------------------------------//
// Test 1: shortopts, no arguments
void tst_shortopts_noargs(rtt_dsxx::UnitTest &ut) {
  std::cout << "\n>>> Begin test 1..." << std::endl;

  // Simulate argc and argv for no arguments case.
  std::array<char const *, 1> my_argv = {"tstXGetopt"};

  // Register options
  rtt_dsxx::XGetopt program_options(static_cast<int>(my_argv.size()), my_argv.data(), "hv");
  int c(0);
  while ((c = program_options()) != -1) {
    switch (c) {
    case 'h':
      ITFAILS;
      break;

    case 'v':
      ITFAILS;
      break;

    default:
      ITFAILS; // no options to parse.
      break;
    }
  }

  FAIL_IF_NOT(program_options.get_option_value() == std::string(""));
  FAIL_IF_NOT(program_options.get_unmatched_arguments().size() == 0);
  PASSMSG("XGetopt with no arguments and short options works as expected.");
  return;
}

//------------------------------------------------------------------------------------------------//
// Test 2: shortopts, version argument.
void tst_shortopts_args_v(rtt_dsxx::UnitTest &ut) {
  std::cout << "\n>>> Begin test 2..." << std::endl;

  // Simulate argc and argv for no arguments case.
  std::array<char const *, 2> my_argv = {"tstXGetopt", "-v"};

  // Register options
  rtt_dsxx::XGetopt program_options(static_cast<int>(my_argv.size()), my_argv.data(), "hv");
  int c(0);
  while ((c = program_options()) != -1) {
    switch (c) {
    case 'h':
      ITFAILS;
      break;

    case 'v': {
      std::ostringstream msg;
      msg << "shortopts_args_v: XGetopt with '-v' argument and short options"
          << " works as expected.";
      PASSMSG(msg.str());
      break;
    }

    default:
      ITFAILS; // no options to parse.
      break;
    }
  }

  FAIL_IF_NOT(program_options.get_option_value() == std::string(""));
  FAIL_IF_NOT(program_options.get_unmatched_arguments().size() == 0);
  return;
}

//------------------------------------------------------------------------------------------------//
// Test 3: shortopts, unknown argument.
void tst_shortopts_args_x(rtt_dsxx::UnitTest &ut) {
  std::cout << "\n>>> Begin test 3..." << std::endl;

  // Simulate argc and argv for no arguments case.
  std::array<char const *, 2> my_argv = {"tstXGetopt", "-x"};

  // Register options
  rtt_dsxx::XGetopt program_options(static_cast<int>(my_argv.size()), my_argv.data(), "hv");
  int c(0);
  while ((c = program_options()) != -1) {
    switch (c) {
    case 'h':
      ITFAILS;
      break;

    case 'v':
      ITFAILS;
      break;

    default:
      ITFAILS; // no options to parse.
      break;
    }
  }

  FAIL_IF_NOT(program_options.get_option_value() == std::string(""));
  FAIL_IF_NOT(program_options.get_unmatched_arguments()[0] == "-x");
  PASSMSG("XGetopt with unknown arguments and short options works as expected.");
  return;
}

//------------------------------------------------------------------------------------------------//
// Test 4: shortopts, 2 args (known and unknown)
void tst_shortopts_args_vx(rtt_dsxx::UnitTest &ut) {
  std::cout << "\n>>> Begin test 4..." << std::endl;

  // Simulate argc and argv for no arguments case.
  std::array<char const *, 3> my_argv = {"tstXGetopt", "-x", "-v"};

  // Register options
  rtt_dsxx::XGetopt program_options(static_cast<int>(my_argv.size()), my_argv.data(), "hv");
  int c(0);
  while ((c = program_options()) != -1) {
    switch (c) {
    case 'h':
      ITFAILS;
      break;

    case 'v': {
      std::ostringstream msg;
      msg << "shortopts_args_vx: XGetopt with '-v' argument and short options"
          << " works as expected.";
      PASSMSG(msg.str());
      break;
    }

    default:
      ITFAILS; // no options to parse.
      break;
    }
  }

  FAIL_IF_NOT(program_options.get_option_value() == std::string(""));
  FAIL_IF_NOT(program_options.get_unmatched_arguments()[0] == "-x");
  std::ostringstream msg;
  msg << "XGetopt with known and unknown arguments and short options works as"
      << " expected.";
  PASSMSG(msg.str());
  return;
}

//------------------------------------------------------------------------------------------------//
// Test 5: shortopts, 2 args (one with value)
void tst_shortopts_args_vc(rtt_dsxx::UnitTest &ut) {
  std::cout << "\n>>> Begin test 5..." << std::endl;

  // Simulate argc and argv for no arguments case.
  std::array<char const *, 4> my_argv = {"tstXGetopt", "-v", "-c", "pink"};

  // Register options
  rtt_dsxx::XGetopt program_options(static_cast<int>(my_argv.size()), my_argv.data(), "hvc:");
  int c(0);
  while ((c = program_options()) != -1) {
    switch (c) {
    case 'c': {
      std::ostringstream msg;
      if (program_options.get_option_value() == "pink") {
        msg << "shortopts_args_vc: Read user specified value correctly.";
        PASSMSG(msg.str());
      } else {
        msg << "shortopts_args_vc: Did not read user specified value correctly";
        FAILMSG(msg.str());
      }
      break;
    }

    case 'h':
      ITFAILS;
      break;

    case 'v': {
      std::ostringstream msg;
      msg << "shortopts_args_vx: XGetopt with '-v' argument and short options"
          << " works as expected.";
      PASSMSG(msg.str());
      break;
    }

    default:
      ITFAILS; // no options to parse.
      break;
    }
  }

  FAIL_IF_NOT(program_options.get_option_value() == std::string("pink"));
  FAIL_IF_NOT(program_options.get_unmatched_arguments().size() == 0);
  return;
}

//------------------------------------------------------------------------------------------------//
// Test 6: shortopts, help message
void tst_shortopts_args_h(rtt_dsxx::UnitTest &ut) {
  std::cout << "\n>>> Begin test 6..." << std::endl;

  // Simulate argc and argv for no arguments case.
  std::array<char const *, 2> my_argv = {"tstXGetopt", "-h"};

  // Register options
  rtt_dsxx::XGetopt program_options(static_cast<int>(my_argv.size()), my_argv.data(), "hvc:");
  int c(0);
  while ((c = program_options()) != -1) {
    switch (c) {
    case 'c':
      ITFAILS;
      break;

    case 'h':
      PASSMSG("shortopts_args_h: help option selected.");
      break;

    case 'v':
      ITFAILS;
      break;

    default:
      ITFAILS; // no options to parse.
      break;
    }
  }

  FAIL_IF_NOT(program_options.get_option_value() == std::string(""));
  FAIL_IF_NOT(program_options.get_unmatched_arguments().size() == 0);

  std::ostringstream helpmsg;
  helpmsg << program_options.display_help("tst_shortopts_args_h");
  FAIL_IF(helpmsg.str().size() == 0);
  std::cout << helpmsg.str() << std::endl;
  bool verbose(true);
  std::map<std::string, unsigned> word_list(rtt_dsxx::get_word_count(helpmsg, verbose));
  FAIL_IF_NOT(word_list[std::string("Options")] == 2);
  FAIL_IF_NOT(word_list[std::string("<value>")] == 1);
  FAIL_IF_NOT(word_list[std::string("-v")] == 1);
  return;
}

//------------------------------------------------------------------------------------------------//
// Test 7: shortopts, double dash
void tst_shortopts_args_doubledash(rtt_dsxx::UnitTest &ut) {
  std::cout << "\n>>> Begin test 7..." << std::endl;

  // Simulate argc and argv for no arguments case.
  std::array<char const *, 4> my_argv = {"tstXGetopt", "-h", "--", "-v"};

  // Register options
  rtt_dsxx::XGetopt program_options(static_cast<int>(my_argv.size()), my_argv.data(), "hvc:");
  int c(0);
  while ((c = program_options()) != -1) {
    switch (c) {
    case 'c':
      ITFAILS;
      break;

    case 'h':
      PASSMSG("shortopts_args_doubledash: matched option -h");
      break;

    case 'v': {
      std::ostringstream msg;
      msg << "shortopts_args_doubledash: version option selected. Should not"
          << " have arrived here.";
      FAILMSG(msg.str());
      break;
    }

    default:
      ITFAILS; // no options to parse.
      break;
    }
  }

  FAIL_IF_NOT(program_options.get_option_value() == std::string(""));
  FAIL_IF_NOT(program_options.get_unmatched_arguments().size() == 0);
  return;
}

//------------------------------------------------------------------------------------------------//
// Test 8: shorotps, missing argument
void tst_shortopts_args_missingarg(rtt_dsxx::UnitTest &ut) {
  std::cout << "\n>>> Begin test 8..." << std::endl;

  // Simulate argc and argv for no arguments case.
  std::array<char const *, 2> my_argv = {"tstXGetopt", "-c"};

  // Register options
  try {
    rtt_dsxx::XGetopt program_options(static_cast<int>(my_argv.size()), my_argv.data(), "hvc:");
    FAILMSG("Insist failed to fire with missing required argument.");
  } catch (rtt_dsxx::assertion & /*error*/) {
    PASSMSG("As expected, Insist fired with missing required argument.");
  }
  return;
}

//------------------------------------------------------------------------------------------------//
// Test 9: shortopts, help message with help strings
void tst_shortopts_args_helpstrings(rtt_dsxx::UnitTest &ut) {
  std::cout << "\n>>> Begin test 9..." << std::endl;

  // Simulate argc and argv for no arguments case.
  std::array<char const *, 2> my_argv = {"tstXGetopt", "-h"};

  // Register options
  std::map<char, std::string> help_strings = {{'c', "create a new file with provided\nfilename."},
                                              {'h', "print this message."},
                                              {'v', "print version information and exit."}};

  rtt_dsxx::XGetopt program_options(static_cast<int>(my_argv.size()), my_argv.data(),
                                    "hvc:", help_strings);
  int c(0);
  while ((c = program_options()) != -1) {
    switch (c) {
    case 'c':
      ITFAILS;
      break;

    case 'h':
      PASSMSG("shortopts_args_h: help option selected.");
      break;

    case 'v':
      ITFAILS;
      break;

    default:
      ITFAILS; // no options to parse.
      break;
    }
  }

  FAIL_IF_NOT(program_options.get_option_value() == std::string(""));
  FAIL_IF_NOT(program_options.get_unmatched_arguments().size() == 0);

  std::ostringstream helpmsg;
  helpmsg << program_options.display_help("tst_shortopts_args_h");
  FAIL_IF(helpmsg.str().size() == 0);
  std::cout << helpmsg.str() << std::endl;
  bool verbose(true);
  std::map<std::string, unsigned> word_list(rtt_dsxx::get_word_count(helpmsg, verbose));
  FAIL_IF_NOT(word_list[std::string("Options")] == 2);
  FAIL_IF_NOT(word_list[std::string("<value>")] == 1);
  FAIL_IF_NOT(word_list[std::string("-v")] == 1);
  FAIL_IF_NOT(word_list[std::string("message")] == 1);
  FAIL_IF_NOT(word_list[std::string("print")] == 2);
  return;
}

//------------------------------------------------------------------------------------------------//
// Test 10: long options, no arguments
void tst_lopts_noargs(rtt_dsxx::UnitTest &ut) {
  std::cout << "\n>>> Begin test 10..." << std::endl;

  // Simulate argc and argv for no arguments case.
  std::array<char const *, 1> my_argv = {"tstXGetopt"};

  // Register options
  rtt_dsxx::XGetopt::csmap long_options = {{'c', "create:"}, {'h', "help"}, {'v', "version"}};

  rtt_dsxx::XGetopt program_options(static_cast<int>(my_argv.size()), my_argv.data(), long_options);
  int c(0);
  while ((c = program_options()) != -1) {
    switch (c) {
    case 'h':
      ITFAILS;
      break;

    case 'v':
      ITFAILS;
      break;

    default:
      ITFAILS; // no options to parse.
      break;
    }
  }

  FAIL_IF_NOT(program_options.get_option_value() == std::string(""));
  FAIL_IF_NOT(program_options.get_unmatched_arguments().size() == 0);
  std::ostringstream msg;
  msg << "tst_lopts_noargs: XGetopt with no arguments and long options works as"
      << " expected.";
  PASSMSG(msg.str());
  return;
}

//------------------------------------------------------------------------------------------------//
// Test 11: longopts, version argument.
void tst_lopts_args_v(rtt_dsxx::UnitTest &ut) {
  std::cout << "\n>>> Begin test 11..." << std::endl;

  // Simulate argc and argv for no arguments case.
  std::array<char const *, 2> my_argv = {"tstXGetopt", "-v"};

  // Register options
  rtt_dsxx::XGetopt::csmap long_options = {{'c', "create:"}, {'h', "help"}, {'v', "version"}};

  rtt_dsxx::XGetopt program_options(static_cast<int>(my_argv.size()), my_argv.data(), long_options);
  int c(0);
  while ((c = program_options()) != -1) {
    switch (c) {
    case 'h':
      ITFAILS;
      break;

    case 'v': {
      std::ostringstream msg;
      msg << "lopts_args_v: XGetopt with '-v' argument and long "
          << "options works as expected.";
      PASSMSG(msg.str());
      break;
    }

    default:
      ITFAILS; // no options to parse.
      break;
    }
  }

  FAIL_IF_NOT(program_options.get_option_value() == std::string(""));
  FAIL_IF_NOT(program_options.get_unmatched_arguments().size() == 0);
  return;
}

//------------------------------------------------------------------------------------------------//
// Test 12: lopts, version argument.
void tst_lopts_args_version(rtt_dsxx::UnitTest &ut) {
  std::cout << "\n>>> Begin test 12..." << std::endl;

  // Simulate argc and argv for no arguments case.
  std::array<char const *, 2> my_argv = {"tstXGetopt", "--version"};

  // Register options
  rtt_dsxx::XGetopt::csmap long_options = {{'c', "create:"}, {'h', "help"}, {'v', "version"}};

  rtt_dsxx::XGetopt program_options(static_cast<int>(my_argv.size()), my_argv.data(), long_options);
  int c(0);
  while ((c = program_options()) != -1) {
    switch (c) {
    case 'h':
      ITFAILS;
      break;

    case 'v': {
      std::ostringstream msg;
      msg << "lopts_args_v: XGetopt with '--version' argument and "
          << "long options works as expected.";
      PASSMSG(msg.str());
      break;
    }

    default:
      ITFAILS; // no options to parse.
      break;
    }
  }

  FAIL_IF_NOT(program_options.get_option_value() == std::string(""));
  FAIL_IF_NOT(program_options.get_unmatched_arguments().size() == 0);
  return;
}

//------------------------------------------------------------------------------------------------//
// Test 13: lopts, unknown argument.
void tst_lopts_args_xray(rtt_dsxx::UnitTest &ut) {
  std::cout << "\n>>> Begin test 13..." << std::endl;

  // Simulate argc and argv for no arguments case.
  std::array<char const *, 3> my_argv = {"tstXGetopt", "--xray", "yellow"};

  // Register options
  rtt_dsxx::XGetopt::csmap long_options = {{'c', "create:"}, {'h', "help"}, {'v', "version"}};

  rtt_dsxx::XGetopt program_options(static_cast<int>(my_argv.size()), my_argv.data(), long_options);
  int c(0);
  while ((c = program_options()) != -1) {
    switch (c) {
    case 'h':
      ITFAILS;
      break;

    case 'v':
      ITFAILS;
      break;

    default:
      ITFAILS; // no options to parse.
      break;
    }
  }

  FAIL_IF_NOT(program_options.get_option_value() == std::string(""));
  FAIL_IF_NOT(program_options.get_unmatched_arguments()[0] == "--xray");
  std::ostringstream msg;
  msg << "XGetopt with unknown arguments and long options works as expected.";
  PASSMSG(msg.str());
  return;
}

//------------------------------------------------------------------------------------------------//
// Test 14: lopts, 2 args (long and short)
void tst_lopts_args_vhelp(rtt_dsxx::UnitTest &ut) {
  std::cout << "\n>>> Begin test 14..." << std::endl;

  // Simulate argc and argv for no arguments case.
  std::array<char const *, 3> my_argv = {"tstXGetopt", "-v", "--help"};

  // Register options
  rtt_dsxx::XGetopt::csmap long_options = {{'c', "create:"}, {'h', "help"}, {'v', "version"}};

  rtt_dsxx::XGetopt program_options(static_cast<int>(my_argv.size()), my_argv.data(), long_options);

  int c(0);
  bool vflag(false);
  bool hflag(false);
  while ((c = program_options()) != -1) {
    switch (c) {
    case 'h':
      hflag = true;
      break;

    case 'v':
      vflag = true;
      break;

    default:
      ITFAILS; // no options to parse.
      break;
    }
  }

  if (hflag && vflag) {
    std::ostringstream msg;
    msg << "lopts_args_vhelp: XGetopt with '-v' and '--help' arguments and long"
        << " options works as expected.";
    PASSMSG(msg.str());
  }
  FAIL_IF_NOT(program_options.get_option_value() == std::string(""));
  FAIL_IF_NOT(program_options.get_unmatched_arguments().size() == 0);

  return;
}

//------------------------------------------------------------------------------------------------//
// Test 15: lopts, 2 args (one with value)
void tst_lopts_args_versioncreate(rtt_dsxx::UnitTest &ut) {
  std::cout << "\n>>> Begin test 15..." << std::endl;

  // Simulate argc and argv for no arguments case.
  std::array<char const *, 4> my_argv = {"tstXGetopt", "--version", "--create", "pink"};

  // Register options
  rtt_dsxx::XGetopt::csmap long_options = {{'c', "create:"}, {'h', "help"}, {'v', "version"}};

  rtt_dsxx::XGetopt program_options(static_cast<int>(my_argv.size()), my_argv.data(), long_options);

  bool cflag(false);
  bool vflag(false);

  int c(0);
  while ((c = program_options()) != -1) {
    switch (c) {
    case 'c':
      cflag = true;
      if (program_options.get_option_value() == "pink") {
        std::ostringstream msg;
        msg << "lopts_args_versioncreate: Read user specified value correctly.";
        PASSMSG(msg.str());
      } else {
        std::ostringstream msg;
        msg << "lopts_args_versioncreate: Did not read user specified value"
            << " correctly.";
        FAILMSG(msg.str());
      }
      break;

    case 'h':
      ITFAILS;
      break;

    case 'v':
      vflag = true;
      break;

    default:
      ITFAILS; // no options to parse.
      break;
    }
  }

  if (cflag && vflag)
    PASSMSG("lopts_args_versioncreate: parsed both c and v options.");
  FAIL_IF_NOT(program_options.get_option_value() == std::string("pink"));
  FAIL_IF_NOT(program_options.get_unmatched_arguments().size() == 0);
  return;
}

//------------------------------------------------------------------------------------------------//
// Test 16: lopts, help message
void tst_lopts_args_help(rtt_dsxx::UnitTest &ut) {
  std::cout << "\n>>> Begin test 16..." << std::endl;

  // Simulate argc and argv for no arguments case.
  std::array<char const *, 2> my_argv = {"tstXGetopt", "--help"};

  // Register options
  rtt_dsxx::XGetopt::csmap long_options = {{'c', "create:"}, {'h', "help"}, {'v', "version"}};

  rtt_dsxx::XGetopt program_options(static_cast<int>(my_argv.size()), my_argv.data(), long_options);

  int c(0);
  while ((c = program_options()) != -1) {
    switch (c) {
    case 'c':
      ITFAILS;
      break;

    case 'h':
      PASSMSG("lopts_args_help: help option selected.");
      break;

    case 'v':
      ITFAILS;
      break;

    default:
      ITFAILS; // no options to parse.
      break;
    }
  }

  FAIL_IF_NOT(program_options.get_option_value() == std::string(""));
  FAIL_IF_NOT(program_options.get_unmatched_arguments().size() == 0);

  std::ostringstream helpmsg;
  helpmsg << program_options.display_help("tst_lopts_args_h");
  FAIL_IF(helpmsg.str().size() == 0);
  std::cout << helpmsg.str() << std::endl;
  bool verbose(true);
  std::map<std::string, unsigned> word_list(rtt_dsxx::get_word_count(helpmsg, verbose));
  FAIL_IF_NOT(word_list[std::string("Options")] == 2);
  FAIL_IF_NOT(word_list[std::string("--help")] == 1);
  FAIL_IF_NOT(word_list[std::string("--create")] == 1);
  return;
}

//------------------------------------------------------------------------------------------------//
// Test 17: lopts, double dash
void tst_lopts_args_doubledash(rtt_dsxx::UnitTest &ut) {
  std::cout << "\n>>> Begin test 17..." << std::endl;

  // Simulate argc and argv for no arguments case.
  std::array<char const *, 4> my_argv = {"tstXGetopt", "--help", "--", "--version"};

  // Register options
  rtt_dsxx::XGetopt::csmap long_options = {{'c', "create:"}, {'h', "help"}, {'v', "version"}};

  rtt_dsxx::XGetopt program_options(static_cast<int>(my_argv.size()), my_argv.data(), long_options);

  int c(0);
  while ((c = program_options()) != -1) {
    switch (c) {
    case 'c':
      ITFAILS;
      break;

    case 'h':
      PASSMSG("lopts_args_doubledash: matched option --help");
      break;

    case 'v': {
      std::ostringstream msg;
      msg << "lopts_args_doubledash: version option selected. Should not have"
          << " arrived here.";
      FAILMSG(msg.str());
      break;
    }

    default:
      ITFAILS; // no options to parse.
      break;
    }
  }

  FAIL_IF_NOT(program_options.get_option_value() == std::string(""));
  FAIL_IF_NOT(program_options.get_unmatched_arguments().size() == 0);
  return;
}

//------------------------------------------------------------------------------------------------//
// Test 18: shorotps, missing argument
void tst_lopts_args_missingarg(rtt_dsxx::UnitTest &ut) {
  std::cout << "\n>>> Begin test 18..." << std::endl;

  // Simulate argc and argv for no arguments case.
  std::array<char const *, 2> my_argv = {"tstXGetopt", "--create"};

  // Register options
  try {
    rtt_dsxx::XGetopt::csmap long_options = {{'c', "create:"}, {'h', "help"}, {'v', "version"}};

    rtt_dsxx::XGetopt program_options(static_cast<int>(my_argv.size()), my_argv.data(),
                                      long_options);

    FAILMSG("Insist failed to fire with missing required argument.");
  } catch (rtt_dsxx::assertion & /*error*/) {
    PASSMSG("As expected, Insist fired with missing required argument.");
  }
  return;
}

//------------------------------------------------------------------------------------------------//
// Test 19: lopts, help message with help strings
void tst_lopts_args_helpstrings(rtt_dsxx::UnitTest &ut) {
  std::cout << "\n>>> Begin test 19..." << std::endl;

  // Simulate argc and argv for no arguments case.
  std::array<char const *, 2> my_argv = {"tstXGetopt", "--help"};

  // Register options
  rtt_dsxx::XGetopt::csmap long_options = {
      {'a', "appendtoit:"}, {'c', "create:"}, {'h', "help"}, {'v', "version"}};

  std::map<char, std::string> help_strings = {
      {'a', "append to an existing file with\nprovided filename."},
      {'c', "create a new file with provided\nfilename."},
      {'h', "print this message."}};
  // test print when missing.
  // help_strings['v'] = "print version information and exit.";

  rtt_dsxx::XGetopt program_options(static_cast<int>(my_argv.size()), my_argv.data(), long_options,
                                    help_strings);
  int c(0);
  while ((c = program_options()) != -1) {
    switch (c) {
    case 'c':
      ITFAILS;
      break;

    case 'h':
      PASSMSG("lopts_args_h: --help option selected.");
      break;

    case 'v':
      ITFAILS;
      break;

    default:
      ITFAILS; // no options to parse.
      break;
    }
  }

  FAIL_IF_NOT(program_options.get_option_value() == std::string(""));
  FAIL_IF_NOT(program_options.get_unmatched_arguments().size() == 0);

  std::ostringstream helpmsg;
  helpmsg << program_options.display_help("tst_lopts_args_helpstrings");
  FAIL_IF(helpmsg.str().size() == 0);
  std::cout << helpmsg.str() << std::endl;
  bool verbose(true);
  std::map<std::string, unsigned> word_list(rtt_dsxx::get_word_count(helpmsg, verbose));
  FAIL_IF_NOT(word_list[std::string("Options")] == 2);
  FAIL_IF_NOT(word_list[std::string("<value>")] == 2);
  FAIL_IF_NOT(word_list[std::string("-v")] == 1);
  FAIL_IF_NOT(word_list[std::string("message")] == 1);
  FAIL_IF_NOT(word_list[std::string("print")] == 1);
  return;
}

//------------------------------------------------------------------------------------------------//
int main(int argc, char *argv[]) {
  rtt_dsxx::ScalarUnitTest ut(argc, argv, rtt_dsxx::release);
  try {
    tst_shortopts_noargs(ut);           // 1
    tst_shortopts_args_v(ut);           // 2
    tst_shortopts_args_x(ut);           // 3
    tst_shortopts_args_vx(ut);          // 4
    tst_shortopts_args_vc(ut);          // 5
    tst_shortopts_args_h(ut);           // 6
    tst_shortopts_args_doubledash(ut);  // 7
    tst_shortopts_args_missingarg(ut);  // 8
    tst_shortopts_args_helpstrings(ut); // 9

    tst_lopts_noargs(ut);             // 10
    tst_lopts_args_v(ut);             // 11
    tst_lopts_args_version(ut);       // 12
    tst_lopts_args_xray(ut);          // 13
    tst_lopts_args_vhelp(ut);         // 14
    tst_lopts_args_versioncreate(ut); // 15
    tst_lopts_args_help(ut);          // 16
    tst_lopts_args_doubledash(ut);    // 17
    tst_lopts_args_missingarg(ut);    // 18
    tst_lopts_args_helpstrings(ut);   // 19
  }
  UT_EPILOG(ut);
}

//------------------------------------------------------------------------------------------------//
// end of tstXGetopt.cc
//------------------------------------------------------------------------------------------------//
