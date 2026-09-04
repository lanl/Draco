/*
Copyright 2016, D. E. Shaw Research.
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are
met:

* Redistributions of source code must retain the above copyright
  notice, this list of conditions, and the following disclaimer.

* Redistributions in binary form must reproduce the above copyright
  notice, this list of conditions, and the following disclaimer in the
  documentation and/or other materials provided with the distribution.

* Neither the name of D. E. Shaw Research nor the names of its
  contributors may be used to endorse or promote products derived from
  this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#ifdef __GNUC__
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wfloat-equal"
#pragma GCC diagnostic ignored "-Wconversion"
#pragma GCC diagnostic ignored "-Wunused-but-set-variable"
#endif

#ifdef __clang__
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wexpansion-to-defined"
#pragma clang diagnostic ignored "-Wreserved-id-macro"
#endif

#ifdef _MSC_FULL_VER
// 4311 - 'type cast': pointer truncation from 'char[54]' to 'long'
// 4202 - nonstandard extension used: non-constant aggregate initializer
// 4127 - conditional expression is constant
// 4302 - 'type cast': truncation from 'const char *' to 'long'
#pragma warning(push)
#pragma warning(disable : 4311 4204 4127 4302)
#endif

#include "rng/config.h"

#include "Random123/conventional/gsl_cbrng.h"
#include "ut_gsl.h"
#include <assert.h>
#include <gsl/gsl_randist.h>
#include <stdio.h>

/* Exercise the GSL_CBRNG macro */

GSL_CBRNG(cbrng, threefry4x64); /* creates gsl_rng_cbrng */

int main(int argc, char **argv) {
  int i;
  gsl_rng *r;
  gsl_rng *rcopy;
  unsigned long x;
  rngRemember(unsigned long save);
  rngRemember(unsigned long saved[5]);
  double sum = 0.;
  (void)argc;
  (void)argv; /* unused */

  /* Silence an unused-parameter warning. */
  (void)argc;

  r = gsl_rng_alloc(gsl_rng_cbrng);
  assert(gsl_rng_min(r) == 0);
  assert(gsl_rng_max(r) == 0xffffffffUL); /* Not necessarily ~0UL */
  assert(gsl_rng_size(r) > 0);

  printf("%s\nulongs from %s in initial state\n", argv[0], gsl_rng_name(r));
  for (i = 0; i < 5; i++) {
    x = gsl_rng_get(r);
    rngRemember(saved[i] = x);
    printf("%d: 0x%lx\n", i, x);
    assert(x != 0);
  }
  printf("uniforms from %s\n", gsl_rng_name(r));
  for (i = 0; i < 5; i++) {
    double z = gsl_rng_uniform(r);
    sum += z;
    printf("%d: %.4g\n", i, z);
  }
  assert(sum < 0.9 * 5 && sum > 0.1 * 5 &&
         (long)"sum must be reasonably close  to 0.5*number of trials");
  rngRemember(save =) gsl_rng_get(r);

  gsl_rng_set(r, 0xdeadbeef); /* set a non-zero seed */
  printf("ulongs from %s after seed\n", gsl_rng_name(r));
  for (i = 0; i < 5; i++) {
    x = gsl_rng_get(r);
    printf("%d: 0x%lx\n", i, x);
    assert(x != 0);
  }
  /* make a copy of the total state */
  rcopy = gsl_rng_alloc(gsl_rng_cbrng);
  gsl_rng_memcpy(rcopy, r);
  printf("uniforms from %s\n", gsl_rng_name(r));
  sum = 0.;
  for (i = 0; i < 5; i++) {
    double x2 = gsl_rng_uniform(r);
    rngRemember(double y = gsl_rng_uniform(rcopy));
    printf("%d: %.4g\n", i, x2);
    sum += x2;
    assert(x2 == y);
  }
  assert(gsl_rng_get(r) != save);
  assert(sum < 0.9 * 5 && sum > 0.1 * 5 &&
         (long)"sum must be reasonably close  to 0.5*number of trials");

  /* gsl_rng_set(*, 0) is supposed to recover the default seed */
  gsl_rng_set(r, 0);
  printf("ulongs from %s after restore to initial\n", gsl_rng_name(r));
  for (i = 0; i < 5; i++) {
    x = gsl_rng_get(r);
    assert(x == saved[i]);
    printf("%d: 0x%lx\n", i, x);
    assert(x != 0);
  }
  printf("uniforms from %s\n", gsl_rng_name(r));
  for (i = 0; i < 5; i++) {
    printf("%d: %.4g\n", i, gsl_rng_uniform(r));
  }

  gsl_rng_free(rcopy);
  gsl_rng_free(r);

  printf("ut_gsl: OK\n");
  return 0;
}

#ifdef _MSC_FULL_VER
#pragma warning(pop)
#endif

#ifdef __clang__
// Restore clang diagnostics to previous state.
#pragma clang diagnostic pop
#endif

#ifdef __GNUC__
// Restore GCC diagnostics to previous state.
#pragma GCC diagnostic pop
#endif
