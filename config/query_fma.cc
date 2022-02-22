//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   config/query_fma.cc
 * \author Kelly Thompson <kgt@lanl.gov>
 * \date   Thursday, Feb 09, 2017, 08:11 am
 * \brief  FMA features test.
 * \note   Copyright (C) 2017-2022 Triad National Security, LLC., All rights reserved.
 *
 * Supporting functions:
 *
 * This code is adopted from
 * https://software.intel.com/en-us/node/405250?language=es&wapkw=avx2+cpuid
 *
 * An alternative approach for Visual Studio can be found at
 * https://docs.microsoft.com/en-us/cpp/intrinsics/cpuid-cpuidex?view=vs-2017
 */
//------------------------------------------------------------------------------------------------//

#include <string>

//------------------------------------------------------------------------------------------------//
// Intel Compiler
//------------------------------------------------------------------------------------------------//

#if defined(__INTEL_COMPILER) && (__INTEL_COMPILER >= 1300)

#include <immintrin.h>

int check_4th_gen_intel_core_features(int const mode = 0) {

  if (mode == 1)
    return _may_i_use_cpu_feature(_FEATURE_AVX2);
  else if (mode == 2)
    return _may_i_use_cpu_feature(_FEATURE_FMA);
  else {
    const int the_4th_gen_features =
        (_FEATURE_AVX2 | _FEATURE_FMA | _FEATURE_BMI | _FEATURE_LZCNT | _FEATURE_MOVBE);
    return _may_i_use_cpu_feature(the_4th_gen_features);
  }
}

//------------------------------------------------------------------------------------------------//
// Non-Intel Compiler
//------------------------------------------------------------------------------------------------//

#else /* non-Intel compiler */

#include <stdint.h>
#if defined(_MSC_VER)
#define MYINT int
#include <intrin.h>
#else
#define MYINT uint32_t
#endif

//------------------------------------------------------------------------------------------------//
void run_cpuid(MYINT eax, MYINT ecx, MYINT *abcd) {
#if defined(_MSC_VER)
  __cpuidex(abcd, eax, ecx);
#else
  MYINT ebx(42), edx(42);
#if defined(__i386__) && defined(__PIC__)
  /* in case of PIC under 32-bit EBX cannot be clobbered */
  __asm__("movl %%ebx, %%edi \n\t cpuid \n\t xchgl %%ebx, %%edi"
          : "=D"(ebx),
#else
  __asm__("cpuid"
          : "+b"(ebx),
#endif
            "+a"(eax), "+c"(ecx), "=d"(edx));
  abcd[0] = eax;
  abcd[1] = ebx;
  abcd[2] = ecx;
  abcd[3] = edx;
#endif
}

//------------------------------------------------------------------------------------------------//
int check_xcr0_ymm() {
  MYINT xcr0;
#if defined(_MSC_VER)
  xcr0 = (MYINT)_xgetbv(0); /* min VS2010 SP1 compiler is required */
#else
  /* named form of xgetbv not supported on OSX, so must use byte form, see:
     https://github.com/asmjit/asmjit/issues/78
   */
  __asm__(".byte 0x0F, 0x01, 0xd0" : "=a"(xcr0) : "c"(0) : "%edx");
#endif
  return ((xcr0 & 6) == 6); /* checking if xmm and ymm state are enabled in XCR0 */
}

//------------------------------------------------------------------------------------------------//
int check_4th_gen_intel_core_features(int const mode = 0) {
  MYINT abcd[4];
  MYINT fma_movbe_osxsave_mask = ((1 << 12) | (1 << 22) | (1 << 27));
  MYINT avx2_bmi12_mask = (1 << 5) | (1 << 3) | (1 << 8);

  bool have_fma(false);
  bool have_avx2(false);

  /* CPUID.(EAX=01H, ECX=0H):ECX.FMA[bit 12]==1   &&
       CPUID.(EAX=01H, ECX=0H):ECX.MOVBE[bit 22]==1 &&
       CPUID.(EAX=01H, ECX=0H):ECX.OSXSAVE[bit 27]==1 */
  run_cpuid(1, 0, abcd);
  if ((abcd[2] & fma_movbe_osxsave_mask) == fma_movbe_osxsave_mask)
    have_fma = true;

  if (mode == 2)
    return have_fma ? 1 : 0;

  /*  CPUID.(EAX=07H, ECX=0H):EBX.AVX2[bit 5]==1  &&
        CPUID.(EAX=07H, ECX=0H):EBX.BMI1[bit 3]==1  &&
        CPUID.(EAX=07H, ECX=0H):EBX.BMI2[bit 8]==1  */
  run_cpuid(7, 0, abcd);
  if ((abcd[1] & avx2_bmi12_mask) == avx2_bmi12_mask)
    have_avx2 = true;

  if (mode == 1)
    return have_avx2 ? 1 : 0;

  if (!check_xcr0_ymm())
    return 0;

  /* CPUID.(EAX=80000001H):ECX.LZCNT[bit 5]==1 */
  run_cpuid(0x80000001, 0, abcd);
  if ((abcd[2] & (1 << 5)) == 0)
    return 0;

  // We have 4 of the 4th-gen cpu features!
  return 1;
}

#endif /* non-Intel compiler */

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Query the CPU for Gen4 features.
 *
 * \return 1 if features are found, else 0.
 *
 * Command line options:
 * -f Only check for FMA
 * -a Only check for AVX2
 */
int main(int argc, char *argv[]) {
  if (argc == 1)
    return check_4th_gen_intel_core_features();
  else {
    for (int i = 1; i < argc; ++i) {
      if (std::string(argv[i]) == std::string("-a"))
        return check_4th_gen_intel_core_features(1);
      if (std::string(argv[i]) == std::string("-f"))
        return check_4th_gen_intel_core_features(2);
    }
  }
}

//------------------------------------------------------------------------------------------------//
// end of query_fma.cc
//------------------------------------------------------------------------------------------------//
