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
void bli_dgemm_rviv_asm_4vx4
    (
      uint64_t             k,
      double*     restrict alpha,
      double*     restrict a,
      double*     restrict b,
      double*     restrict beta,
      double*     restrict c, uint64_t rs_c, uint64_t cs_c
    );

void bli_dgemm_rviv_4vx4
     (
       dim_t               m,
       dim_t               n,
       dim_t               k,
       double*    restrict alpha,
       double*    restrict a,
       double*    restrict b,
       double*    restrict beta,
       double*    restrict c, inc_t rs_c0, inc_t cs_c0,
       auxinfo_t*          data,
       cntx_t*             cntx
     )
{
    // Use local copy in case dim_t has a different size than expected in the assembly kernel
    uint64_t _k     = k;
    uint64_t rs_c   = rs_c0;
    uint64_t cs_c   = cs_c0;

    // Extract vector-length dependent mr, nr that are fixed at configure time.
    const inc_t mr = bli_cntx_get_blksz_def_dt( BLIS_DOUBLE, BLIS_MR, cntx );
    const inc_t nr = 4;

    GEMM_UKR_SETUP_CT( d, mr, nr, false );

    bli_dgemm_rviv_asm_4vx4(_k, alpha, a, b, beta, c, rs_c, cs_c);

    GEMM_UKR_FLUSH_CT( d );
}