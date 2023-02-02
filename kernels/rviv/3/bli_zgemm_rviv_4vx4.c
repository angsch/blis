/*

   BLIS
   An object-based framework for developing high-performance BLAS-like
   libraries.

   Copyright (C) 2023, The University of Texas at Austin

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

#include "bli_rviv_utils.h"

void bli_zgemm_rviv_asm_4vx4
    (
      intptr_t               k,
      dcomplex*     restrict alpha,
      dcomplex*     restrict a,
      dcomplex*     restrict b,
      dcomplex*     restrict beta,
      dcomplex*     restrict c, intptr_t rs_c, intptr_t cs_c
    );


void bli_zgemm_rviv_4vx4
     (
       dim_t               m,
       dim_t               n,
       dim_t               k,
       dcomplex*  restrict alpha,
       dcomplex*  restrict a,
       dcomplex*  restrict b,
       dcomplex*  restrict beta,
       dcomplex*  restrict c, inc_t rs_c, inc_t cs_c,
       auxinfo_t*          data,
       cntx_t*             cntx
     )
{
    // The assembly kernels always take native machine-sized integer arguments.
    // dim_t and inc_t are normally defined as being machine-sized. If larger, assert.
    bli_static_assert( sizeof(dim_t) <= sizeof(intptr_t) &&
                       sizeof(inc_t) <= sizeof(intptr_t) );

    // Extract vector-length dependent mr, nr that are fixed at configure time.
    const inc_t mr = bli_cntx_get_blksz_def_dt( BLIS_DCOMPLEX, BLIS_MR, cntx );
    const inc_t nr = 4;

    GEMM_UKR_SETUP_CT( z, mr, nr, false );

    bli_zgemm_rviv_asm_4vx4( k, alpha, a, b, beta, c,
                             rs_c * get_vlenb(), cs_c * sizeof(dcomplex) );

    GEMM_UKR_FLUSH_CT( z );
}
