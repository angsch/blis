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

	.text
	.align      2
	.global     REALNAME

// void REALNAME(intptr_t k, void* alpha, void* a, void* b,
//               void* beta, void* c, intptr_t rs_c, intptr_t cs_c)
//
// register arguments:
// a0   k
// a1   alpha
// a2   a
// a3   b
// a4   beta
// a5   c
// a6   rs_c
// a7   cs_c
//

#define REALSIZE (DATASIZE/2)

#define loop_counter a0

#define A00_ptr   a2
#define A10_ptr   t0
#define A01_ptr   t1
#define A11_ptr   t2

#define B_row_ptr a3

#define C00_ptr   a5
#define C01_ptr   t3
#define C02_ptr   t4
#define C03_ptr   t5
#define C04_ptr   t0
#define C05_ptr   t1
#define C10_ptr   s1
#define C11_ptr   s2
#define C12_ptr   s3
#define C13_ptr   s4
#define C14_ptr   s5
#define C15_ptr   s6

#define tmp       t6

#define ALPHA_re  fa0
#define ALPHA_im  fa1
#define BETA_re   fa2
#define BETA_im   fa3

#define B00_re    fa4
#define B00_im    fa5
#define B01_re    fa6
#define B01_im    fa7
#define B02_re    fa0
#define B02_im    fa1
#define B03_re    fa2
#define B03_im    fa3
#define B04_re    ft8
#define B04_im    ft9
#define B05_re    ft10
#define B05_im    ft11

#define B10_re    ft0
#define B10_im    ft1
#define B11_re    ft2
#define B11_im    ft3
#define B12_re    ft4
#define B12_im    ft5
#define B13_re    ft6
#define B13_im    ft7
#define B14_re    fa4
#define B14_im    fa5
#define B15_re    fa6
#define B15_im    fa7

#define fzero     ft8

#define A00_re    v24
#define A00_im    v25
#define A10_re    v26
#define A10_im    v27
#define A01_re    v28
#define A01_im    v29
#define A11_re    v30
#define A11_im    v31

// AB, C and tmp are aliased for readability.
#define C0_re     v24
#define C0_im     v25
#define C1_re     v26
#define C1_im     v27
#define C2_re     v28
#define C2_im     v29
#define C3_re     v30
#define C3_im     v31
#define C4_re     v0
#define C4_im     v1
#define C5_re     v2
#define C5_im     v3

#define AB00_re   v0
#define AB00_im   v1
#define AB01_re   v2
#define AB01_im   v3
#define AB02_re   v4
#define AB02_im   v5
#define AB03_re   v6
#define AB03_im   v7
#define AB04_re   v16
#define AB04_im   v17
#define AB05_re   v18
#define AB05_im   v19
#define AB10_re   v8
#define AB10_im   v9
#define AB11_re   v10
#define AB11_im   v11
#define AB12_re   v12
#define AB12_im   v13
#define AB13_re   v14
#define AB13_im   v15
#define AB14_re   v20
#define AB14_im   v21
#define AB15_re   v22
#define AB15_im   v23

#define tmp0_re   v24
#define tmp0_im   v25
#define tmp1_re   v26
#define tmp1_im   v27
#define tmp2_re   v28
#define tmp2_im   v29
#define tmp3_re   v30
#define tmp3_im   v31

#define rs_c  a6
#define cs_c  a7

REALNAME:
	#include "rviv_save_registers.h"

	vsetvli s0, zero, VTYPE, m1, ta, ma
	csrr s0, vlenb
	slli s0, s0, 1

	// Zero-initialize accumulators
	vxor.vv AB00_re, AB00_re, AB00_re
	vxor.vv AB00_im, AB00_im, AB00_im
	vxor.vv AB01_re, AB01_re, AB01_re
	vxor.vv AB01_im, AB01_im, AB01_im
	vxor.vv AB02_re, AB02_re, AB02_re
	vxor.vv AB02_im, AB02_im, AB02_im
	vxor.vv AB03_re, AB03_re, AB03_re
	vxor.vv AB03_im, AB03_im, AB03_im
	vxor.vv AB04_re, AB04_re, AB04_re
	vxor.vv AB04_im, AB04_im, AB04_im
	vxor.vv AB05_re, AB05_re, AB05_re
	vxor.vv AB05_im, AB05_im, AB05_im
	vxor.vv AB10_re, AB10_re, AB10_re
	vxor.vv AB10_im, AB10_im, AB10_im
	vxor.vv AB11_re, AB11_re, AB11_re
	vxor.vv AB11_im, AB11_im, AB11_im
	vxor.vv AB12_re, AB12_re, AB12_re
	vxor.vv AB12_im, AB12_im, AB12_im
	vxor.vv AB13_re, AB13_re, AB13_re
	vxor.vv AB13_im, AB13_im, AB13_im
	vxor.vv AB14_re, AB14_re, AB14_re
	vxor.vv AB14_im, AB14_im, AB14_im
	vxor.vv AB15_re, AB15_re, AB15_re
	vxor.vv AB15_im, AB15_im, AB15_im

	// Handle k == 0
	beqz loop_counter, MULTIPLYBETA

	add A10_ptr, A00_ptr, s0
	slli s0, s0, 1      // length of a column of A in bytes
	add A01_ptr, A00_ptr, s0
	add A11_ptr, A10_ptr, s0

	li tmp, 3
	ble loop_counter, tmp, TAIL_UNROLL_2

	// Preload A and B
	// Load and deinterleave A(:,l)
	VLE A00_re, (A00_ptr)
	VLE A10_re, (A10_ptr)

	// Load B(l,0:5)
	FLOAD B00_re, 0*REALSIZE(B_row_ptr)
	FLOAD B00_im, 1*REALSIZE(B_row_ptr)
	FLOAD B01_re, 2*REALSIZE(B_row_ptr)
	FLOAD B01_im, 3*REALSIZE(B_row_ptr)
	FLOAD B02_re, 4*REALSIZE(B_row_ptr)
	FLOAD B02_im, 5*REALSIZE(B_row_ptr)
	FLOAD B03_re, 6*REALSIZE(B_row_ptr)
	FLOAD B03_im, 7*REALSIZE(B_row_ptr)
	FLOAD B04_re, 8*REALSIZE(B_row_ptr)
	FLOAD B04_im, 9*REALSIZE(B_row_ptr)
	FLOAD B05_re, 10*REALSIZE(B_row_ptr)
	FLOAD B05_im, 11*REALSIZE(B_row_ptr)

	// Load and deinterleave A(:,l+1)
	VLE A01_re, (A01_ptr)
	VLE A11_re, (A11_ptr)

LOOP_UNROLL_4: // loop_counter >= 4
	addi loop_counter, loop_counter, -4

	vfmacc.vf  AB00_re, B00_re, A00_re   // AB(:,0) += A(:,l) * B(l,0)
	vfnmsac.vf AB00_re, B00_im, A00_im
	vfmacc.vf  AB00_im, B00_re, A00_im
	vfmacc.vf  AB00_im, B00_im, A00_re
	vfmacc.vf  AB10_re, B00_re, A10_re
	vfnmsac.vf AB10_re, B00_im, A10_im
	vfmacc.vf  AB10_im, B00_re, A10_im
	vfmacc.vf  AB10_im, B00_im, A10_re

	vfmacc.vf  AB01_re, B01_re, A00_re   // AB(:,1) += A(:,l) * B(l,1)
	vfnmsac.vf AB01_re, B01_im, A00_im
	vfmacc.vf  AB01_im, B01_re, A00_im
	vfmacc.vf  AB01_im, B01_im, A00_re
	vfmacc.vf  AB11_re, B01_re, A10_re
	vfnmsac.vf AB11_re, B01_im, A10_im
	vfmacc.vf  AB11_im, B01_re, A10_im
	vfmacc.vf  AB11_im, B01_im, A10_re

	// Point to A(:,l+2), A(:,l+3)
	add A00_ptr, A01_ptr, s0
	add A10_ptr, A11_ptr, s0
	add A01_ptr, A00_ptr, s0
	add A11_ptr, A10_ptr, s0

	// Load B(l+1,0:3)
	FLOAD B10_re, 12*REALSIZE(B_row_ptr)
	FLOAD B10_im, 13*REALSIZE(B_row_ptr)
	FLOAD B11_re, 14*REALSIZE(B_row_ptr)
	FLOAD B11_im, 15*REALSIZE(B_row_ptr)
	FLOAD B12_re, 16*REALSIZE(B_row_ptr)
	FLOAD B12_im, 17*REALSIZE(B_row_ptr)
	FLOAD B13_re, 18*REALSIZE(B_row_ptr)
	FLOAD B13_im, 19*REALSIZE(B_row_ptr)

	vfmacc.vf  AB02_re, B02_re, A00_re   // AB(:,2) += A(:,l) * B(l,2)
	vfnmsac.vf AB02_re, B02_im, A00_im
	vfmacc.vf  AB02_im, B02_re, A00_im
	vfmacc.vf  AB02_im, B02_im, A00_re
	vfmacc.vf  AB12_re, B02_re, A10_re
	vfnmsac.vf AB12_re, B02_im, A10_im
	vfmacc.vf  AB12_im, B02_re, A10_im
	vfmacc.vf  AB12_im, B02_im, A10_re

	vfmacc.vf  AB03_re, B03_re, A00_re   // AB(:,3) += A(:,l) * B(l,3)
	vfnmsac.vf AB03_re, B03_im, A00_im
	vfmacc.vf  AB03_im, B03_re, A00_im
	vfmacc.vf  AB03_im, B03_im, A00_re
	vfmacc.vf  AB13_re, B03_re, A10_re
	vfnmsac.vf AB13_re, B03_im, A10_im
	vfmacc.vf  AB13_im, B03_re, A10_im
	vfmacc.vf  AB13_im, B03_im, A10_re

	// Load B(l+1,4:5)
	FLOAD B14_re, 20*REALSIZE(B_row_ptr)
	FLOAD B14_im, 21*REALSIZE(B_row_ptr)
	FLOAD B15_re, 22*REALSIZE(B_row_ptr)
	FLOAD B15_im, 23*REALSIZE(B_row_ptr)
	addi B_row_ptr, B_row_ptr, 24*REALSIZE

	vfmacc.vf  AB00_re, B10_re, A01_re   // AB(:,0) += A(:,l+1) * B(l+1,0)
	vfnmsac.vf AB00_re, B10_im, A01_im
	vfmacc.vf  AB00_im, B10_re, A01_im
	vfmacc.vf  AB00_im, B10_im, A01_re
	vfmacc.vf  AB10_re, B10_re, A11_re
	vfnmsac.vf AB10_re, B10_im, A11_im
	vfmacc.vf  AB10_im, B10_re, A11_im
	vfmacc.vf  AB10_im, B10_im, A11_re

	vfmacc.vf  AB04_re, B04_re, A00_re   // AB(:,4) += A(:,l) * B(l,4)
	vfnmsac.vf AB04_re, B04_im, A00_im
	vfmacc.vf  AB04_im, B04_re, A00_im
	vfmacc.vf  AB04_im, B04_im, A00_re
	vfmacc.vf  AB14_re, B04_re, A10_re
	vfnmsac.vf AB14_re, B04_im, A10_im
	vfmacc.vf  AB14_im, B04_re, A10_im
	vfmacc.vf  AB14_im, B04_im, A10_re

	vfmacc.vf  AB04_re, B14_re, A01_re   // AB(:,4) += A(:,l+1) * B(l+1,4)
	vfnmsac.vf AB04_re, B14_im, A01_im
	vfmacc.vf  AB04_im, B14_re, A01_im
	vfmacc.vf  AB04_im, B14_im, A01_re
	vfmacc.vf  AB14_re, B14_re, A11_re
	vfnmsac.vf AB14_re, B14_im, A11_im
	vfmacc.vf  AB14_im, B14_re, A11_im
	vfmacc.vf  AB14_im, B14_im, A11_re

	vfmacc.vf  AB05_re, B05_re, A00_re   // AB(:,5) += A(:,l) * B(l,5)
	vfnmsac.vf AB05_re, B05_im, A00_im
	vfmacc.vf  AB05_im, B05_re, A00_im
	vfmacc.vf  AB05_im, B05_im, A00_re
	vfmacc.vf  AB15_re, B05_re, A10_re
	vfnmsac.vf AB15_re, B05_im, A10_im
	vfmacc.vf  AB15_im, B05_re, A10_im
	vfmacc.vf  AB15_im, B05_im, A10_re

	vfmacc.vf  AB05_re, B15_re, A01_re   // AB(:,5) += A(:,l+1) * B(l+1,5)
	vfnmsac.vf AB05_re, B15_im, A01_im
	vfmacc.vf  AB05_im, B15_re, A01_im
	vfmacc.vf  AB05_im, B15_im, A01_re
	vfmacc.vf  AB15_re, B15_re, A11_re
	vfnmsac.vf AB15_re, B15_im, A11_im
	vfmacc.vf  AB15_im, B15_re, A11_im
	vfmacc.vf  AB15_im, B15_im, A11_re

	// Load and deinterleave A(:,l+2)
	VLE A00_re, (A00_ptr)
	VLE A10_re, (A10_ptr)

	// Load B(l+2, 0:3)
	FLOAD B00_re, 0*REALSIZE(B_row_ptr)
	FLOAD B00_im, 1*REALSIZE(B_row_ptr)
	FLOAD B01_re, 2*REALSIZE(B_row_ptr)
	FLOAD B01_im, 3*REALSIZE(B_row_ptr)
	FLOAD B02_re, 4*REALSIZE(B_row_ptr)
	FLOAD B02_im, 5*REALSIZE(B_row_ptr)
	FLOAD B03_re, 6*REALSIZE(B_row_ptr)
	FLOAD B03_im, 7*REALSIZE(B_row_ptr)

	vfmacc.vf  AB01_re, B11_re, A01_re   // AB(:,1) += A(:,l+1) * B(l+1,1)
	vfnmsac.vf AB01_re, B11_im, A01_im
	vfmacc.vf  AB01_im, B11_re, A01_im
	vfmacc.vf  AB01_im, B11_im, A01_re
	vfmacc.vf  AB11_re, B11_re, A11_re
	vfnmsac.vf AB11_re, B11_im, A11_im
	vfmacc.vf  AB11_im, B11_re, A11_im
	vfmacc.vf  AB11_im, B11_im, A11_re

	vfmacc.vf  AB02_re, B12_re, A01_re   // AB(:,2) += A(:,l+1) * B(l+1,2)
	vfnmsac.vf AB02_re, B12_im, A01_im
	vfmacc.vf  AB02_im, B12_re, A01_im
	vfmacc.vf  AB02_im, B12_im, A01_re
	vfmacc.vf  AB12_re, B12_re, A11_re
	vfnmsac.vf AB12_re, B12_im, A11_im
	vfmacc.vf  AB12_im, B12_re, A11_im
	vfmacc.vf  AB12_im, B12_im, A11_re

	// Load B(l+2, 4:5)
	FLOAD B04_re, 8*REALSIZE(B_row_ptr)
	FLOAD B04_im, 9*REALSIZE(B_row_ptr)
	FLOAD B05_re, 10*REALSIZE(B_row_ptr)
	FLOAD B05_im, 11*REALSIZE(B_row_ptr)

	vfmacc.vf  AB03_re, B13_re, A01_re   // AB(:,3) += A(:,l+1) * B(l+1,3)
	vfnmsac.vf AB03_re, B13_im, A01_im
	vfmacc.vf  AB03_im, B13_re, A01_im
	vfmacc.vf  AB03_im, B13_im, A01_re
	vfmacc.vf  AB13_re, B13_re, A11_re
	vfnmsac.vf AB13_re, B13_im, A11_im
	vfmacc.vf  AB13_im, B13_re, A11_im
	vfmacc.vf  AB13_im, B13_im, A11_re

	// Load and deinterleave A(:,l+3)
	VLE A01_re, (A01_ptr)
	VLE A11_re, (A11_ptr)

	// Point to A(:,l+4), A(:,l+5)
	add A00_ptr, A01_ptr, s0
	add A10_ptr, A11_ptr, s0
	add A01_ptr, A00_ptr, s0
	add A11_ptr, A10_ptr, s0

	// Load B(l+3, 0:3)
	FLOAD B10_re, 12*REALSIZE(B_row_ptr)
	FLOAD B10_im, 13*REALSIZE(B_row_ptr)
	FLOAD B11_re, 14*REALSIZE(B_row_ptr)
	FLOAD B11_im, 15*REALSIZE(B_row_ptr)
	FLOAD B12_re, 16*REALSIZE(B_row_ptr)
	FLOAD B12_im, 17*REALSIZE(B_row_ptr)
	FLOAD B13_re, 18*REALSIZE(B_row_ptr)
	FLOAD B13_im, 19*REALSIZE(B_row_ptr)

	vfmacc.vf  AB00_re, B00_re, A00_re   // AB(:,0) += A(:,l+2) * B(l+2,0)
	vfnmsac.vf AB00_re, B00_im, A00_im
	vfmacc.vf  AB00_im, B00_re, A00_im
	vfmacc.vf  AB00_im, B00_im, A00_re
	vfmacc.vf  AB10_re, B00_re, A10_re
	vfnmsac.vf AB10_re, B00_im, A10_im
	vfmacc.vf  AB10_im, B00_re, A10_im
	vfmacc.vf  AB10_im, B00_im, A10_re

	vfmacc.vf  AB01_re, B01_re, A00_re   // AB(:,1) += A(:,l+2) * B(l+2,1)
	vfnmsac.vf AB01_re, B01_im, A00_im
	vfmacc.vf  AB01_im, B01_re, A00_im
	vfmacc.vf  AB01_im, B01_im, A00_re
	vfmacc.vf  AB11_re, B01_re, A10_re
	vfnmsac.vf AB11_re, B01_im, A10_im
	vfmacc.vf  AB11_im, B01_re, A10_im
	vfmacc.vf  AB11_im, B01_im, A10_re

	vfmacc.vf  AB00_re, B10_re, A01_re   // AB(:,0) += A(:,l+3) * B(l+3,0)
	vfnmsac.vf AB00_re, B10_im, A01_im
	vfmacc.vf  AB00_im, B10_re, A01_im
	vfmacc.vf  AB00_im, B10_im, A01_re
	vfmacc.vf  AB10_re, B10_re, A11_re
	vfnmsac.vf AB10_re, B10_im, A11_im
	vfmacc.vf  AB10_im, B10_re, A11_im
	vfmacc.vf  AB10_im, B10_im, A11_re

	vfmacc.vf  AB01_re, B11_re, A01_re   // AB(:,1) += A(:,l+3) * B(l+3,1)
	vfnmsac.vf AB01_re, B11_im, A01_im
	vfmacc.vf  AB01_im, B11_re, A01_im
	vfmacc.vf  AB01_im, B11_im, A01_re
	vfmacc.vf  AB11_re, B11_re, A11_re
	vfnmsac.vf AB11_re, B11_im, A11_im
	vfmacc.vf  AB11_im, B11_re, A11_im
	vfmacc.vf  AB11_im, B11_im, A11_re

	vfmacc.vf  AB02_re, B02_re, A00_re   // AB(:,2) += A(:,l+2) * B(l+2,2)
	vfnmsac.vf AB02_re, B02_im, A00_im
	vfmacc.vf  AB02_im, B02_re, A00_im
	vfmacc.vf  AB02_im, B02_im, A00_re
	vfmacc.vf  AB12_re, B02_re, A10_re
	vfnmsac.vf AB12_re, B02_im, A10_im
	vfmacc.vf  AB12_im, B02_re, A10_im
	vfmacc.vf  AB12_im, B02_im, A10_re

	// Load B(l+3, 4:5)
	FLOAD B14_re, 20*REALSIZE(B_row_ptr)
	FLOAD B14_im, 21*REALSIZE(B_row_ptr)
	FLOAD B15_re, 22*REALSIZE(B_row_ptr)
	FLOAD B15_im, 23*REALSIZE(B_row_ptr)
	addi B_row_ptr, B_row_ptr, 24*REALSIZE

	vfmacc.vf  AB02_re, B12_re, A01_re   // AB(:,2) += A(:,l+3) * B(l+3,2)
	vfnmsac.vf AB02_re, B12_im, A01_im
	vfmacc.vf  AB02_im, B12_re, A01_im
	vfmacc.vf  AB02_im, B12_im, A01_re
	vfmacc.vf  AB12_re, B12_re, A11_re
	vfnmsac.vf AB12_re, B12_im, A11_im
	vfmacc.vf  AB12_im, B12_re, A11_im
	vfmacc.vf  AB12_im, B12_im, A11_re

	vfmacc.vf  AB03_re, B03_re, A00_re   // AB(:,3) += A(:,l+2) * B(l+2,3)
	vfnmsac.vf AB03_re, B03_im, A00_im
	vfmacc.vf  AB03_im, B03_re, A00_im
	vfmacc.vf  AB03_im, B03_im, A00_re
	vfmacc.vf  AB13_re, B03_re, A10_re
	vfnmsac.vf AB13_re, B03_im, A10_im
	vfmacc.vf  AB13_im, B03_re, A10_im
	vfmacc.vf  AB13_im, B03_im, A10_re

	vfmacc.vf  AB03_re, B13_re, A01_re   // AB(:,3) += A(:,l+3) * B(l+3,3)
	vfnmsac.vf AB03_re, B13_im, A01_im
	vfmacc.vf  AB03_im, B13_re, A01_im
	vfmacc.vf  AB03_im, B13_im, A01_re
	vfmacc.vf  AB13_re, B13_re, A11_re
	vfnmsac.vf AB13_re, B13_im, A11_im
	vfmacc.vf  AB13_im, B13_re, A11_im
	vfmacc.vf  AB13_im, B13_im, A11_re

	vfmacc.vf  AB04_re, B04_re, A00_re   // AB(:,4) += A(:,l+2) * B(l+2,4)
	vfnmsac.vf AB04_re, B04_im, A00_im
	vfmacc.vf  AB04_im, B04_re, A00_im
	vfmacc.vf  AB04_im, B04_im, A00_re
	vfmacc.vf  AB14_re, B04_re, A10_re
	vfnmsac.vf AB14_re, B04_im, A10_im
	vfmacc.vf  AB14_im, B04_re, A10_im
	vfmacc.vf  AB14_im, B04_im, A10_re

	vfmacc.vf  AB04_re, B14_re, A01_re   // AB(:,4) += A(:,l+3) * B(l+3,4)
	vfnmsac.vf AB04_re, B14_im, A01_im
	vfmacc.vf  AB04_im, B14_re, A01_im
	vfmacc.vf  AB04_im, B14_im, A01_re
	vfmacc.vf  AB14_re, B14_re, A11_re
	vfnmsac.vf AB14_re, B14_im, A11_im
	vfmacc.vf  AB14_im, B14_re, A11_im
	vfmacc.vf  AB14_im, B14_im, A11_re

	vfmacc.vf  AB05_re, B05_re, A00_re   // AB(:,5) += A(:,l+2) * B(l+2,5)
	vfnmsac.vf AB05_re, B05_im, A00_im
	vfmacc.vf  AB05_im, B05_re, A00_im
	vfmacc.vf  AB05_im, B05_im, A00_re
	vfmacc.vf  AB15_re, B05_re, A10_re
	vfnmsac.vf AB15_re, B05_im, A10_im
	vfmacc.vf  AB15_im, B05_re, A10_im
	vfmacc.vf  AB15_im, B05_im, A10_re

	vfmacc.vf  AB05_re, B15_re, A01_re   // AB(:,5) += A(:,l+3) * B(l+3,5)
	vfnmsac.vf AB05_re, B15_im, A01_im
	vfmacc.vf  AB05_im, B15_re, A01_im
	vfmacc.vf  AB05_im, B15_im, A01_re
	vfmacc.vf  AB15_re, B15_re, A11_re
	vfnmsac.vf AB15_re, B15_im, A11_im
	vfmacc.vf  AB15_im, B15_re, A11_im
	vfmacc.vf  AB15_im, B15_im, A11_re

	li tmp, 3
	ble loop_counter, tmp, TAIL_UNROLL_2

	// Load A and B for the next iteration
	VLE A00_re, (A00_ptr)
	VLE A10_re, (A10_ptr)
	VLE A01_re, (A01_ptr)
	VLE A11_re, (A11_ptr)

	FLOAD B00_re, 0*REALSIZE(B_row_ptr)
	FLOAD B00_im, 1*REALSIZE(B_row_ptr)
	FLOAD B01_re, 2*REALSIZE(B_row_ptr)
	FLOAD B01_im, 3*REALSIZE(B_row_ptr)
	FLOAD B02_re, 4*REALSIZE(B_row_ptr)
	FLOAD B02_im, 5*REALSIZE(B_row_ptr)
	FLOAD B03_re, 6*REALSIZE(B_row_ptr)
	FLOAD B03_im, 7*REALSIZE(B_row_ptr)
	FLOAD B04_re, 8*REALSIZE(B_row_ptr)
	FLOAD B04_im, 9*REALSIZE(B_row_ptr)
	FLOAD B05_re, 10*REALSIZE(B_row_ptr)
	FLOAD B05_im, 11*REALSIZE(B_row_ptr)

	j LOOP_UNROLL_4

TAIL_UNROLL_2: // loop_counter <= 3
	li tmp, 1
	ble loop_counter, tmp, TAIL_UNROLL_1

	addi loop_counter, loop_counter, -2

	// Load and deinterleave A(:,l)
	VLE A00_re, (A00_ptr)
	VLE A10_re, (A10_ptr)

	// Load B(l, 0:3)
	FLOAD B00_re, 0*REALSIZE(B_row_ptr)
	FLOAD B00_im, 1*REALSIZE(B_row_ptr)
	FLOAD B01_re, 2*REALSIZE(B_row_ptr)
	FLOAD B01_im, 3*REALSIZE(B_row_ptr)
	FLOAD B02_re, 4*REALSIZE(B_row_ptr)
	FLOAD B02_im, 5*REALSIZE(B_row_ptr)
	FLOAD B03_re, 6*REALSIZE(B_row_ptr)
	FLOAD B03_im, 7*REALSIZE(B_row_ptr)

	vfmacc.vf  AB00_re, B00_re, A00_re   // AB(:,0) += A(:,l) * B(l,0)
	vfnmsac.vf AB00_re, B00_im, A00_im
	vfmacc.vf  AB00_im, B00_re, A00_im
	vfmacc.vf  AB00_im, B00_im, A00_re
	vfmacc.vf  AB10_re, B00_re, A10_re
	vfnmsac.vf AB10_re, B00_im, A10_im
	vfmacc.vf  AB10_im, B00_re, A10_im
	vfmacc.vf  AB10_im, B00_im, A10_re

	vfmacc.vf  AB01_re, B01_re, A00_re   // AB(:,1) += A(:,l) * B(l,1)
	vfnmsac.vf AB01_re, B01_im, A00_im
	vfmacc.vf  AB01_im, B01_re, A00_im
	vfmacc.vf  AB01_im, B01_im, A00_re
	vfmacc.vf  AB11_re, B01_re, A10_re
	vfnmsac.vf AB11_re, B01_im, A10_im
	vfmacc.vf  AB11_im, B01_re, A10_im
	vfmacc.vf  AB11_im, B01_im, A10_re

	// Load and deinterleave A(:,l+1)
	VLE A01_re, (A01_ptr)
	VLE A11_re, (A11_ptr)

	// Load B(l+1, 0:3)
	FLOAD B10_re, 12*REALSIZE(B_row_ptr)
	FLOAD B10_im, 13*REALSIZE(B_row_ptr)
	FLOAD B11_re, 14*REALSIZE(B_row_ptr)
	FLOAD B11_im, 15*REALSIZE(B_row_ptr)
	FLOAD B12_re, 16*REALSIZE(B_row_ptr)
	FLOAD B12_im, 17*REALSIZE(B_row_ptr)
	FLOAD B13_re, 18*REALSIZE(B_row_ptr)
	FLOAD B13_im, 19*REALSIZE(B_row_ptr)

	vfmacc.vf  AB00_re, B10_re, A01_re   // AB(:,0) += A(:,l+1) * B(l+1,0)
	vfnmsac.vf AB00_re, B10_im, A01_im
	vfmacc.vf  AB00_im, B10_re, A01_im
	vfmacc.vf  AB00_im, B10_im, A01_re
	vfmacc.vf  AB10_re, B10_re, A11_re
	vfnmsac.vf AB10_re, B10_im, A11_im
	vfmacc.vf  AB10_im, B10_re, A11_im
	vfmacc.vf  AB10_im, B10_im, A11_re

	vfmacc.vf  AB01_re, B11_re, A01_re   // AB(:,1) += A(:,l+1) * B(l+1,1)
	vfnmsac.vf AB01_re, B11_im, A01_im
	vfmacc.vf  AB01_im, B11_re, A01_im
	vfmacc.vf  AB01_im, B11_im, A01_re
	vfmacc.vf  AB11_re, B11_re, A11_re
	vfnmsac.vf AB11_re, B11_im, A11_im
	vfmacc.vf  AB11_im, B11_re, A11_im
	vfmacc.vf  AB11_im, B11_im, A11_re

	vfmacc.vf  AB02_re, B02_re, A00_re   // AB(:,2) += A(:,l) * B(l,2)
	vfnmsac.vf AB02_re, B02_im, A00_im
	vfmacc.vf  AB02_im, B02_re, A00_im
	vfmacc.vf  AB02_im, B02_im, A00_re
	vfmacc.vf  AB12_re, B02_re, A10_re
	vfnmsac.vf AB12_re, B02_im, A10_im
	vfmacc.vf  AB12_im, B02_re, A10_im
	vfmacc.vf  AB12_im, B02_im, A10_re

	vfmacc.vf  AB03_re, B03_re, A00_re   // AB(:,3) += A(:,l) * B(l,3)
	vfnmsac.vf AB03_re, B03_im, A00_im
	vfmacc.vf  AB03_im, B03_re, A00_im
	vfmacc.vf  AB03_im, B03_im, A00_re
	vfmacc.vf  AB13_re, B03_re, A10_re
	vfnmsac.vf AB13_re, B03_im, A10_im
	vfmacc.vf  AB13_im, B03_re, A10_im
	vfmacc.vf  AB13_im, B03_im, A10_re

	vfmacc.vf  AB02_re, B12_re, A01_re   // AB(:,2) += A(:,l+1) * B(l+1,2)
	vfnmsac.vf AB02_re, B12_im, A01_im
	vfmacc.vf  AB02_im, B12_re, A01_im
	vfmacc.vf  AB02_im, B12_im, A01_re
	vfmacc.vf  AB12_re, B12_re, A11_re
	vfnmsac.vf AB12_re, B12_im, A11_im
	vfmacc.vf  AB12_im, B12_re, A11_im
	vfmacc.vf  AB12_im, B12_im, A11_re

	vfmacc.vf  AB03_re, B13_re, A01_re   // AB(:,3) += A(:,l+1) * B(l+1,3)
	vfnmsac.vf AB03_re, B13_im, A01_im
	vfmacc.vf  AB03_im, B13_re, A01_im
	vfmacc.vf  AB03_im, B13_im, A01_re
	vfmacc.vf  AB13_re, B13_re, A11_re
	vfnmsac.vf AB13_re, B13_im, A11_im
	vfmacc.vf  AB13_im, B13_re, A11_im
	vfmacc.vf  AB13_im, B13_im, A11_re

	// Load B(l, 4:5)
	FLOAD B04_re, 8*REALSIZE(B_row_ptr)
	FLOAD B04_im, 9*REALSIZE(B_row_ptr)
	FLOAD B05_re, 10*REALSIZE(B_row_ptr)
	FLOAD B05_im, 11*REALSIZE(B_row_ptr)

	vfmacc.vf  AB04_re, B04_re, A00_re   // AB(:,4) += A(:,l) * B(l,4)
	vfnmsac.vf AB04_re, B04_im, A00_im
	vfmacc.vf  AB04_im, B04_re, A00_im
	vfmacc.vf  AB04_im, B04_im, A00_re
	vfmacc.vf  AB14_re, B04_re, A10_re
	vfnmsac.vf AB14_re, B04_im, A10_im
	vfmacc.vf  AB14_im, B04_re, A10_im
	vfmacc.vf  AB14_im, B04_im, A10_re

	vfmacc.vf  AB05_re, B05_re, A00_re   // AB(:,5) += A(:,l) * B(l,5)
	vfnmsac.vf AB05_re, B05_im, A00_im
	vfmacc.vf  AB05_im, B05_re, A00_im
	vfmacc.vf  AB05_im, B05_im, A00_re
	vfmacc.vf  AB15_re, B05_re, A10_re
	vfnmsac.vf AB15_re, B05_im, A10_im
	vfmacc.vf  AB15_im, B05_re, A10_im
	vfmacc.vf  AB15_im, B05_im, A10_re

	// Load B(l+1, 4:5)
	FLOAD B14_re, 20*REALSIZE(B_row_ptr)
	FLOAD B14_im, 21*REALSIZE(B_row_ptr)
	FLOAD B15_re, 22*REALSIZE(B_row_ptr)
	FLOAD B15_im, 23*REALSIZE(B_row_ptr)

	vfmacc.vf  AB04_re, B14_re, A01_re   // AB(:,4) += A(:,l+1) * B(l+1,4)
	vfnmsac.vf AB04_re, B14_im, A01_im
	vfmacc.vf  AB04_im, B14_re, A01_im
	vfmacc.vf  AB04_im, B14_im, A01_re
	vfmacc.vf  AB14_re, B14_re, A11_re
	vfnmsac.vf AB14_re, B14_im, A11_im
	vfmacc.vf  AB14_im, B14_re, A11_im
	vfmacc.vf  AB14_im, B14_im, A11_re

	vfmacc.vf  AB05_re, B15_re, A01_re   // AB(:,5) += A(:,l+1) * B(l+1,5)
	vfnmsac.vf AB05_re, B15_im, A01_im
	vfmacc.vf  AB05_im, B15_re, A01_im
	vfmacc.vf  AB05_im, B15_im, A01_re
	vfmacc.vf  AB15_re, B15_re, A11_re
	vfnmsac.vf AB15_re, B15_im, A11_im
	vfmacc.vf  AB15_im, B15_re, A11_im
	vfmacc.vf  AB15_im, B15_im, A11_re

	beqz loop_counter, MULTIPLYALPHA

	// Advance pointers
	add A00_ptr, A01_ptr, s0
	add A10_ptr, A11_ptr, s0
	addi B_row_ptr, B_row_ptr, 24*REALSIZE

TAIL_UNROLL_1: // loop_counter <= 1
	beqz loop_counter, MULTIPLYALPHA

	// Load and deinterleave A(:,l)
	VLE A00_re, (A00_ptr)
	VLE A10_re, (A10_ptr)

	// Load B(l,0:3)
	FLOAD B00_re, 0*REALSIZE(B_row_ptr)
	FLOAD B00_im, 1*REALSIZE(B_row_ptr)
	FLOAD B01_re, 2*REALSIZE(B_row_ptr)
	FLOAD B01_im, 3*REALSIZE(B_row_ptr)
	FLOAD B02_re, 4*REALSIZE(B_row_ptr)
	FLOAD B02_im, 5*REALSIZE(B_row_ptr)
	FLOAD B03_re, 6*REALSIZE(B_row_ptr)
	FLOAD B03_im, 7*REALSIZE(B_row_ptr)

	vfmacc.vf  AB00_re, B00_re, A00_re   // AB(:,0) += A(:,l) * B(l,0)
	vfnmsac.vf AB00_re, B00_im, A00_im
	vfmacc.vf  AB00_im, B00_re, A00_im
	vfmacc.vf  AB00_im, B00_im, A00_re
	vfmacc.vf  AB10_re, B00_re, A10_re
	vfnmsac.vf AB10_re, B00_im, A10_im
	vfmacc.vf  AB10_im, B00_re, A10_im
	vfmacc.vf  AB10_im, B00_im, A10_re

	vfmacc.vf  AB01_re, B01_re, A00_re   // AB(:,1) += A(:,l) * B(l,1)
	vfnmsac.vf AB01_re, B01_im, A00_im
	vfmacc.vf  AB01_im, B01_re, A00_im
	vfmacc.vf  AB01_im, B01_im, A00_re
	vfmacc.vf  AB11_re, B01_re, A10_re
	vfnmsac.vf AB11_re, B01_im, A10_im
	vfmacc.vf  AB11_im, B01_re, A10_im
	vfmacc.vf  AB11_im, B01_im, A10_re

	// Load B(l,4:5)
	FLOAD B04_re, 8*REALSIZE(B_row_ptr)
	FLOAD B04_im, 9*REALSIZE(B_row_ptr)
	FLOAD B05_re, 10*REALSIZE(B_row_ptr)
	FLOAD B05_im, 11*REALSIZE(B_row_ptr)

	vfmacc.vf  AB02_re, B02_re, A00_re   // AB(:,2) += A(:,l) * B(l,2)
	vfnmsac.vf AB02_re, B02_im, A00_im
	vfmacc.vf  AB02_im, B02_re, A00_im
	vfmacc.vf  AB02_im, B02_im, A00_re
	vfmacc.vf  AB12_re, B02_re, A10_re
	vfnmsac.vf AB12_re, B02_im, A10_im
	vfmacc.vf  AB12_im, B02_re, A10_im
	vfmacc.vf  AB12_im, B02_im, A10_re

	vfmacc.vf  AB03_re, B03_re, A00_re   // AB(:,3) += A(:,l) * B(l,3)
	vfnmsac.vf AB03_re, B03_im, A00_im
	vfmacc.vf  AB03_im, B03_re, A00_im
	vfmacc.vf  AB03_im, B03_im, A00_re
	vfmacc.vf  AB13_re, B03_re, A10_re
	vfnmsac.vf AB13_re, B03_im, A10_im
	vfmacc.vf  AB13_im, B03_re, A10_im
	vfmacc.vf  AB13_im, B03_im, A10_re

	vfmacc.vf  AB04_re, B04_re, A00_re   // AB(:,4) += A(:,l) * B(l,4)
	vfnmsac.vf AB04_re, B04_im, A00_im
	vfmacc.vf  AB04_im, B04_re, A00_im
	vfmacc.vf  AB04_im, B04_im, A00_re
	vfmacc.vf  AB14_re, B04_re, A10_re
	vfnmsac.vf AB14_re, B04_im, A10_im
	vfmacc.vf  AB14_im, B04_re, A10_im
	vfmacc.vf  AB14_im, B04_im, A10_re

	vfmacc.vf  AB05_re, B05_re, A00_re   // AB(:,5) += A(:,l) * B(l,5)
	vfnmsac.vf AB05_re, B05_im, A00_im
	vfmacc.vf  AB05_im, B05_re, A00_im
	vfmacc.vf  AB05_im, B05_im, A00_re
	vfmacc.vf  AB15_re, B05_re, A10_re
	vfnmsac.vf AB15_re, B05_im, A10_im
	vfmacc.vf  AB15_im, B05_re, A10_im
	vfmacc.vf  AB15_im, B05_im, A10_re

MULTIPLYALPHA:
	FLOAD ALPHA_re, 0*REALSIZE(a1)
	FLOAD ALPHA_im, 1*REALSIZE(a1)
	FZERO(fzero)

	FEQ tmp, ALPHA_im, fzero
	bne tmp, zero, ALPHAREAL

	// [AB00, ..., AB03] * alpha
	vfmul.vf  tmp0_re, AB00_im, ALPHA_im
	vfmul.vf  tmp0_im, AB00_re, ALPHA_im
	vfmul.vf  tmp1_re, AB01_im, ALPHA_im
	vfmul.vf  tmp1_im, AB01_re, ALPHA_im
	vfmul.vf  tmp2_re, AB02_im, ALPHA_im
	vfmul.vf  tmp2_im, AB02_re, ALPHA_im
	vfmul.vf  tmp3_re, AB03_im, ALPHA_im
	vfmul.vf  tmp3_im, AB03_re, ALPHA_im
	vfmsub.vf AB00_re, ALPHA_re, tmp0_re
	vfmsub.vf AB01_re, ALPHA_re, tmp1_re
	vfmsub.vf AB02_re, ALPHA_re, tmp2_re
	vfmsub.vf AB03_re, ALPHA_re, tmp3_re
	vfmadd.vf AB00_im, ALPHA_re, tmp0_im
	vfmadd.vf AB01_im, ALPHA_re, tmp1_im
	vfmadd.vf AB02_im, ALPHA_re, tmp2_im
	vfmadd.vf AB03_im, ALPHA_re, tmp3_im

	// [AB10, ..., AB13] * alpha
	vfmul.vf  tmp0_re, AB10_im, ALPHA_im
	vfmul.vf  tmp0_im, AB10_re, ALPHA_im
	vfmul.vf  tmp1_re, AB11_im, ALPHA_im
	vfmul.vf  tmp1_im, AB11_re, ALPHA_im
	vfmul.vf  tmp2_re, AB12_im, ALPHA_im
	vfmul.vf  tmp2_im, AB12_re, ALPHA_im
	vfmul.vf  tmp3_re, AB13_im, ALPHA_im
	vfmul.vf  tmp3_im, AB13_re, ALPHA_im
	vfmsub.vf AB10_re, ALPHA_re, tmp0_re
	vfmsub.vf AB11_re, ALPHA_re, tmp1_re
	vfmsub.vf AB12_re, ALPHA_re, tmp2_re
	vfmsub.vf AB13_re, ALPHA_re, tmp3_re
	vfmadd.vf AB10_im, ALPHA_re, tmp0_im
	vfmadd.vf AB11_im, ALPHA_re, tmp1_im
	vfmadd.vf AB12_im, ALPHA_re, tmp2_im
	vfmadd.vf AB13_im, ALPHA_re, tmp3_im

	// [AB04, AB05] * alpha
	vfmul.vf  tmp0_re, AB04_im, ALPHA_im
	vfmul.vf  tmp0_im, AB04_re, ALPHA_im
	vfmsub.vf AB04_re, ALPHA_re, tmp0_re
	vfmadd.vf AB04_im, ALPHA_re, tmp0_im
	vfmul.vf  tmp1_re, AB05_im, ALPHA_im
	vfmul.vf  tmp1_im, AB05_re, ALPHA_im
	vfmsub.vf AB05_re, ALPHA_re, tmp1_re
	vfmadd.vf AB05_im, ALPHA_re, tmp1_im

	// [AB14, AB15] * alpha
	vfmul.vf  tmp2_re, AB14_im, ALPHA_im
	vfmul.vf  tmp2_im, AB14_re, ALPHA_im
	vfmsub.vf AB14_re, ALPHA_re, tmp2_re
	vfmadd.vf AB14_im, ALPHA_re, tmp2_im
	vfmul.vf  tmp3_re, AB15_im, ALPHA_im
	vfmul.vf  tmp3_im, AB15_re, ALPHA_im
	vfmsub.vf AB15_re, ALPHA_re, tmp3_re
	vfmadd.vf AB15_im, ALPHA_re, tmp3_im
	j MULTIPLYBETA

ALPHAREAL:
	vfmul.vf AB00_re, AB00_re, ALPHA_re
	vfmul.vf AB00_im, AB00_im, ALPHA_re
	vfmul.vf AB01_re, AB01_re, ALPHA_re
	vfmul.vf AB01_im, AB01_im, ALPHA_re
	vfmul.vf AB02_re, AB02_re, ALPHA_re
	vfmul.vf AB02_im, AB02_im, ALPHA_re
	vfmul.vf AB03_re, AB03_re, ALPHA_re
	vfmul.vf AB03_im, AB03_im, ALPHA_re
	vfmul.vf AB04_re, AB04_re, ALPHA_re
	vfmul.vf AB04_im, AB04_im, ALPHA_re
	vfmul.vf AB05_re, AB05_re, ALPHA_re
	vfmul.vf AB05_im, AB05_im, ALPHA_re

	vfmul.vf AB10_re, AB10_re, ALPHA_re
	vfmul.vf AB10_im, AB10_im, ALPHA_re
	vfmul.vf AB11_re, AB11_re, ALPHA_re
	vfmul.vf AB11_im, AB11_im, ALPHA_re
	vfmul.vf AB12_re, AB12_re, ALPHA_re
	vfmul.vf AB12_im, AB12_im, ALPHA_re
	vfmul.vf AB13_re, AB13_re, ALPHA_re
	vfmul.vf AB13_im, AB13_im, ALPHA_re
	vfmul.vf AB14_re, AB14_re, ALPHA_re
	vfmul.vf AB14_im, AB14_im, ALPHA_re
	vfmul.vf AB15_re, AB15_re, ALPHA_re
	vfmul.vf AB15_im, AB15_im, ALPHA_re

MULTIPLYBETA:
	FZERO(fzero)
	FLOAD BETA_re,  0*REALSIZE(a4)
	FLOAD BETA_im,  1*REALSIZE(a4)

	// Set up pointers
	add C01_ptr, C00_ptr, cs_c
	add C02_ptr, C01_ptr, cs_c
	add C03_ptr, C02_ptr, cs_c
	add C04_ptr, C03_ptr, cs_c
	add C05_ptr, C04_ptr, cs_c
	add C10_ptr, C00_ptr, rs_c
	add C11_ptr, C01_ptr, rs_c
	add C12_ptr, C02_ptr, rs_c
	add C13_ptr, C03_ptr, rs_c
	add C14_ptr, C04_ptr, rs_c
	add C15_ptr, C05_ptr, rs_c

	FEQ tmp, BETA_im, fzero
	bne tmp, zero, BETAREAL

	// Load and deinterleave C(0:VLEN-1, 0:1)
	VLE C0_re, (C00_ptr)
	VLE C1_re, (C01_ptr)

	// Load and deinterleave C(0:VLEN-1, 2:3)
	VLE C2_re, (C02_ptr)
	VLE C3_re, (C03_ptr)

	// C(0:VLEN-1,0:1) * beta + AB(0:VLEN-1,0:1)
	vfmacc.vf   AB00_re, BETA_re, C0_re
	vfnmsac.vf  AB00_re, BETA_im, C0_im
	vfmacc.vf   AB00_im, BETA_re, C0_im
	vfmacc.vf   AB00_im, BETA_im, C0_re
	VSE AB00_re, (C00_ptr)

	vfmacc.vf   AB01_re, BETA_re, C1_re
	vfnmsac.vf  AB01_re, BETA_im, C1_im
	vfmacc.vf   AB01_im, BETA_re, C1_im
	vfmacc.vf   AB01_im, BETA_im, C1_re
	VSE AB01_re, (C01_ptr)

	// Load and deinterleave C(0:VLEN-1, 4:5)
	VLE C4_re, (C04_ptr)
	VLE C5_re, (C05_ptr)

	// C(0:VLEN-1,2:3) * beta + AB(0:VLEN-1,2:3)
	vfmacc.vf   AB02_re, BETA_re, C2_re
	vfnmsac.vf  AB02_re, BETA_im, C2_im
	vfmacc.vf   AB02_im, BETA_re, C2_im
	vfmacc.vf   AB02_im, BETA_im, C2_re
	VSE AB02_re, (C02_ptr)

	vfmacc.vf   AB03_re, BETA_re, C3_re
	vfnmsac.vf  AB03_re, BETA_im, C3_im
	vfmacc.vf   AB03_im, BETA_re, C3_im
	vfmacc.vf   AB03_im, BETA_im, C3_re
	VSE AB03_re, (C03_ptr)

	// C(0:VLEN-1,4:5) * beta + AB(0:VLEN-1,4:5)
	vfmacc.vf   AB04_re, BETA_re, C4_re
	vfnmsac.vf  AB04_re, BETA_im, C4_im
	vfmacc.vf   AB04_im, BETA_re, C4_im
	vfmacc.vf   AB04_im, BETA_im, C4_re
	VSE AB04_re, (C04_ptr)

	vfmacc.vf   AB05_re, BETA_re, C5_re
	vfnmsac.vf  AB05_re, BETA_im, C5_im
	vfmacc.vf   AB05_im, BETA_re, C5_im
	vfmacc.vf   AB05_im, BETA_im, C5_re
	VSE AB05_re, (C05_ptr)

	// Load and deinterleave C(VLEN:2*VLEN-1, 0:1)
	VLE C0_re, (C10_ptr)
	VLE C1_re, (C11_ptr)

	// Load and deinterleave C(VLEN:2*VLEN-1, 2:3)
	VLE C2_re, (C12_ptr)
	VLE C3_re, (C13_ptr)

	// C(VLEN:2*VLEN-1,0:1) * beta + AB(VLEN:2*VLEN-1,0:1)
	vfmacc.vf   AB10_re, BETA_re, C0_re
	vfnmsac.vf  AB10_re, BETA_im, C0_im
	vfmacc.vf   AB10_im, BETA_re, C0_im
	vfmacc.vf   AB10_im, BETA_im, C0_re
	VSE AB10_re, (C10_ptr)

	vfmacc.vf   AB11_re, BETA_re, C1_re
	vfnmsac.vf  AB11_re, BETA_im, C1_im
	vfmacc.vf   AB11_im, BETA_re, C1_im
	vfmacc.vf   AB11_im, BETA_im, C1_re
	VSE AB11_re, (C11_ptr)

	// Load and deinterleave C(VLEN:2*VLEN-1, 4:5)
	VLE C4_re, (C14_ptr)
	VLE C5_re, (C15_ptr)

	// C(VLEN:2*VLEN-1,2:3) * beta + AB(VLEN:2*VLEN-1,2:3)
	vfmacc.vf   AB12_re, BETA_re, C2_re
	vfnmsac.vf  AB12_re, BETA_im, C2_im
	vfmacc.vf   AB12_im, BETA_re, C2_im
	vfmacc.vf   AB12_im, BETA_im, C2_re
	VSE AB12_re, (C12_ptr)

	vfmacc.vf   AB13_re, BETA_re, C3_re
	vfnmsac.vf  AB13_re, BETA_im, C3_im
	vfmacc.vf   AB13_im, BETA_re, C3_im
	vfmacc.vf   AB13_im, BETA_im, C3_re
	VSE AB13_re, (C13_ptr)

	// C(VLEN:2*VLEN-1,4:5) * beta + AB(VLEN:2*VLEN-1,4:5)
	vfmacc.vf   AB14_re, BETA_re, C4_re
	vfnmsac.vf  AB14_re, BETA_im, C4_im
	vfmacc.vf   AB14_im, BETA_re, C4_im
	vfmacc.vf   AB14_im, BETA_im, C4_re
	VSE AB14_re, (C14_ptr)

	vfmacc.vf   AB15_re, BETA_re, C5_re
	vfnmsac.vf  AB15_re, BETA_im, C5_im
	vfmacc.vf   AB15_im, BETA_re, C5_im
	vfmacc.vf   AB15_im, BETA_im, C5_re
	VSE AB15_re, (C15_ptr)

	j END

BETAREAL:
	FEQ tmp, BETA_re, fzero
	bne tmp, zero, BETAZERO

	// Load and deinterleave C(0:VLEN-1, 0:3)
	VLE C0_re, (C00_ptr)
	VLE C1_re, (C01_ptr)
	VLE C2_re, (C02_ptr)
	VLE C3_re, (C03_ptr)

	// C(0:VLEN-1,0:3) * beta + AB(0:VLEN-1,0:3)
	vfmacc.vf   AB00_re, BETA_re, C0_re
	vfmacc.vf   AB00_im, BETA_re, C0_im
	vfmacc.vf   AB01_re, BETA_re, C1_re
	vfmacc.vf   AB01_im, BETA_re, C1_im

	vfmacc.vf   AB02_re, BETA_re, C2_re
	vfmacc.vf   AB02_im, BETA_re, C2_im
	vfmacc.vf   AB03_re, BETA_re, C3_re
	vfmacc.vf   AB03_im, BETA_re, C3_im

	VSE AB00_re, (C00_ptr)
	VSE AB01_re, (C01_ptr)
	VSE AB02_re, (C02_ptr)
	VSE AB03_re, (C03_ptr)

	// Load and deinterleave C(0:VLEN-1, 4:5)
	VLE C4_re, (C04_ptr)
	VLE C5_re, (C05_ptr)

	// C(0:VLEN-1,4:5) * beta + AB(0:VLEN-1,4:5)
	vfmacc.vf   AB04_re, BETA_re, C4_re
	vfmacc.vf   AB04_im, BETA_re, C4_im
	vfmacc.vf   AB05_re, BETA_re, C5_re
	vfmacc.vf   AB05_im, BETA_re, C5_im

	VSE AB04_re, (C04_ptr)
	VSE AB05_re, (C05_ptr)

	// Load and deinterleave C(VLEN:2*VLEN-1, 0:3)
	VLE C0_re, (C10_ptr)
	VLE C1_re, (C11_ptr)
	VLE C2_re, (C12_ptr)
	VLE C3_re, (C13_ptr)

	// C(VLEN:2*VLEN-1,0:3) * beta + AB(VLEN:2*VLEN-1,0:3)
	vfmacc.vf   AB10_re, BETA_re, C0_re
	vfmacc.vf   AB10_im, BETA_re, C0_im
	vfmacc.vf   AB11_re, BETA_re, C1_re
	vfmacc.vf   AB11_im, BETA_re, C1_im

	vfmacc.vf   AB12_re, BETA_re, C2_re
	vfmacc.vf   AB12_im, BETA_re, C2_im
	vfmacc.vf   AB13_re, BETA_re, C3_re
	vfmacc.vf   AB13_im, BETA_re, C3_im

	VSE AB10_re, (C10_ptr)
	VSE AB11_re, (C11_ptr)
	VSE AB12_re, (C12_ptr)
	VSE AB13_re, (C13_ptr)

	// Load and deinterleave C(VLEN:2*VLEN-1, 4:5)
	VLE C4_re, (C14_ptr)
	VLE C5_re, (C15_ptr)

	// C(VLEN:2*VLEN-1,4:5) * beta + AB(VLEN:2*VLEN-1,4:5)
	vfmacc.vf   AB14_re, BETA_re, C4_re
	vfmacc.vf   AB14_im, BETA_re, C4_im
	vfmacc.vf   AB15_re, BETA_re, C5_re
	vfmacc.vf   AB15_im, BETA_re, C5_im

	VSE AB14_re, (C14_ptr)
	VSE AB15_re, (C15_ptr)

	j END

BETAZERO:
	VSE AB00_re, (C00_ptr)
	VSE AB01_re, (C01_ptr)
	VSE AB02_re, (C02_ptr)
	VSE AB03_re, (C03_ptr)
	VSE AB04_re, (C04_ptr)
	VSE AB05_re, (C05_ptr)

	VSE AB10_re, (C10_ptr)
	VSE AB11_re, (C11_ptr)
	VSE AB12_re, (C12_ptr)
	VSE AB13_re, (C13_ptr)
	VSE AB14_re, (C14_ptr)
	VSE AB15_re, (C15_ptr)

END:
	#include "rviv_restore_registers.h"
	ret
