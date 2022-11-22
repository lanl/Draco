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
// Check our AES implementation against the example in FIPS-197

#include "rng/config.h"

#ifdef _MSC_FULL_VER
// - 4521: Engines have multiple copy constructors, quite legal C++, disable
//         MSVC complaint.
// - 4244: possible loss of data when converting between int types.
// - 4204: nonstandard extension used - non-constant aggregate initializer
// - 4127: conditional expression is constant
// - 4100: unreferenced formal parameter
#pragma warning(push)
#pragma warning(disable : 4521 4244 4127 4100)
#endif

#ifdef __clang__
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wreserved-id-macro"
#endif

#include "ut_aes.hh"

#include <Random123/ReinterpretCtr.hpp>
#if R123_USE_AES_OPENSSL
#include <openssl/aes.h>
#endif
#include <cassert>
#include <cstdio>
#include <cstring>
#include <iostream>
#include <string>

using namespace std;
using namespace r123;

#if !R123_USE_SSE
int main(int, char **) {
  std::cout << "No SSE support.  This test is not compiled.  OK\n";
  return 0;
}
#else

#include "util_m128.h"

int main(int /*unused*/, char ** /*unused*/) {
#if R123_USE_AES_NI || R123_USE_AES_OPENSSL
  r123array1xm128i IN, K;

  K.v[0].m = m128i_from_charbuf("0001020304050607 08090a0b0c0d0e0f");
  IN.v[0].m = m128i_from_charbuf("0011223344556677 8899aabbccddeeff");
#endif
  // From FIPS-197, this is the official "right answer"
  r123array1xm128i right_answer; // NOLINT [hicpp-member-init]
                                 // nvhpc can't process empty brace initializer here.
  right_answer[0] = m128i_from_charbuf("69c4 e0d8 6a7b 0430 d8cd b780 70b4 c55a");
  (void)
      right_answer; /* don't complain about an unused variable if neither NI nor OPENSSL are enabled. */
#if R123_USE_AES_NI
  if (haveAESNI()) {
    AESNI1xm128i::key_type xk(K);
    AESNI1xm128i bx;
    AESNI1xm128i::ctr_type x = bx(IN, xk);

    assert(x == right_answer);
    cout << "IN: " << m128i_to_string(IN[0]) << "\n";
    cout << "K : " << m128i_to_string(K[0]) << "\n";
    cout << "AES:" << m128i_to_string(x[0]) << "\n";
    cout << "Hooray!  AESNI1xm128i(IN, K) matches the published test vector!\n";
  } else {
    cout << "The AES-NI instructions are not available on this hardware.  "
            "Skipping AES-NI tests\n";
  }
#else
  cout << "The AES-NI Bijections are not compiled into this binary.  Skipping "
          "AES-NI tests\n";
#endif

  // And let's do it with AESOpenSSL.  But since AESOpenSSL has its own
  // format for keys and counters we make a union for the key types and
  // use ReinterpretCtr to wrap a union around the counter types.
#if R123_USE_AES_OPENSSL
#if R123_USE_AES_NI
  typedef AESNI1xm128i::ctr_type nictype;
#else
  typedef r123array1xm128i nictype;
#endif
  AESOpenSSL16x8::ukey_type ouk;
  _mm_storeu_si128((__m128i *)&ouk.v[0], K.v[0].m);
  AESOpenSSL16x8::key_type okey(ouk);
  ReinterpretCtr<nictype, AESOpenSSL16x8> osslb;
  assert(osslb(IN, okey) == right_answer);
  cout << "Hooray!  AESOpenSSL16x8(IN, K) matches the published test vector!\n";
#else
  cout << "The OpenSSL AES implementation is not linked with this binary.  "
          "Skipping the AESOpenSSL16x8\n";
#endif // R123_USE_AES_OPENSSL

  cout << "ut_aes: all OK" << endl;
  return 0;
}

#endif

#ifdef __clang__
// Restore clang diagnostics to previous state.
#pragma clang diagnostic pop
#endif

#ifdef _MSC_FULL_VER
#pragma warning(pop)
#endif
