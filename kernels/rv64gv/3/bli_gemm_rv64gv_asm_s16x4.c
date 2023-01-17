/*

   BLIS
   An object-based framework for developing high-performance BLAS-like
   libraries.

   Copyright (C) 2014, The University of Texas at Austin

   Redistribution and use in source and binary forms, with or without
   modification, are permitted provided that the following conditions are
   met:
	- Redistributions of source code must retain the above copyright
	  notice, this list of conditions and the following disclaimer.
	- Redistributions in binary form must reproduce the above copyright
	  notice, this list of conditions and the following disclaimer in the
	  documentation and/or other materials provided with the distribution.
	- Neither the name(s) of the copyright holder(s) nor the names of its
	  contributors may be used to endorse or promote products derived
	  from this software without specific prior written permission.

   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
   "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
   LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
   A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
   HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
   SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
   LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
   DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
   THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
   (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
   OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


*/
#include "blis.h"

extern
void bli_sgemm_rv64gv_ker_16x4
	 (
	   uint64_t			   k,
	   float*	  restrict alpha,
	   float*	  restrict a,
	   float*	  restrict b,
	   float*	  restrict beta,
	   float*	  restrict c, uint64_t rs_c, uint64_t cs_c
	 );


static uint64_t num_fp32_per_vector() {
	uint64_t velem = 0;
	__asm__ volatile
	(
	" li a5, 256						 \n\t"
	" vsetvli a4, a5, e32, m1, ta, ma	 \n\t"
	" mv %[velem], a4					 \n\t"
	: [velem] "=r" (velem)
	:
	: "a4", "a5"
	);
	return velem;
}

void bli_sgemm_rv64gv_asm_16x4
	 (
	   dim_t			   m,
	   dim_t			   n,
	   dim_t			   k,
	   float*	  restrict alpha,
	   float*	  restrict a,
	   float*	  restrict b,
	   float*	  restrict beta,
	   float*	  restrict c, inc_t rs_c0, inc_t cs_c0,
	   auxinfo_t*		   data,
	   cntx_t*			   cntx
	 )
{
	// Use local copy in case dim_t has a different size than expected in the assembly kernel
	uint64_t _k		= k;
	uint64_t rs_c	= rs_c0;
	uint64_t cs_c	= cs_c0;

	GEMM_UKR_SETUP_CT( s, 16, 4, false );

	// The assembly kernel assumes that the vector length is 128 bits.
	// For now, fall back to generic implementation
	// if VLEN != 128.
	uint64_t velem = num_fp32_per_vector();

	if (velem == 4) {
		bli_sgemm_rv64gv_ker_16x4(_k, alpha, a, b, beta, c, rs_c, cs_c);
	}
	else {
		float ab[16 * 4];
		for (dim_t i = 0; i < 16 * 4; i++) {
			ab[i] = 0;
		}

		if (*alpha) {
			for (dim_t l = 0; l < k; ++l) {
				for (dim_t j = 0; j < 4; j++) {
					for (dim_t i = 0; i < 16; ++i) {
						float bj = b[j];
						float ai = a[i];
						ab[i + j * 16] += ai * bj;
					}
				}
				a += 16;
				b += 4;
			}
			for (dim_t i = 0; i < 16 * 4; i++) {
				ab[i] = ab[i] * (*alpha);
			}
		}

		if (*beta) {
			for (dim_t j = 0; j < n; j++) {
				for (dim_t i = 0; i < m; i++) {
					c[i * rs_c + j * cs_c] = c[i * rs_c + j * cs_c] * (*beta) + ab[i + j * 16];
				}
			}
		}
		else {
			for (dim_t j = 0; j < n; j++) {
				for (dim_t i = 0; i < m; i++) {
					c[i * rs_c + j * cs_c] = ab[i + j * 16];
				}
			}
		}
	}
	GEMM_UKR_FLUSH_CT( s );
}
