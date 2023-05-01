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


#define loop_counter a0

#define A00_ptr   a2
#define A10_ptr   t0
#define A01_ptr   t1
#define A11_ptr   t2
#define A02_ptr   s5
#define A12_ptr   s6
#define A03_ptr   s7
#define A13_ptr   t6

#define B_row_ptr a3

#define C0_ptr    a5
#define C1_ptr    t3
#define C2_ptr    t4
#define C3_ptr    t5
#define C4_ptr    s1
#define C5_ptr    s2
#define C6_ptr    s3
#define C7_ptr    s4
#define C8_ptr    s5
#define C9_ptr    s6

#define tmp    t6

#define ALPHA  fa1
#define BETA   fa2

#define B00    fa0
#define B10    fa1
#define B01    fa2
#define B11    fa3
#define B02    fa4
#define B12    fa5
#define B03    fa6
#define B13    fa7
#define B04    ft0
#define B14    ft1
#define B05    ft2
#define B15    ft3
#define B06    ft4
#define B16    ft5
#define B07    ft6
#define B17    ft7
#define B08    ft8
#define B18    ft9
#define B09    ft10
#define B19    ft11

// aliases of B00, ..., B19 for easier readability
#define B20    fa0
#define B30    fa1
#define B21    fa2
#define B31    fa3
#define B22    fa4
#define B32    fa5
#define B23    fa6
#define B33    fa7
#define B24    ft0
#define B34    ft1
#define B25    ft2
#define B35    ft3
#define B26    ft4
#define B36    ft5
#define B27    ft6
#define B37    ft7
#define B28    ft8
#define B38    ft9
#define B29    ft10
#define B39    ft11

#define fzero  ft8

#define A00    v20
#define A10    v21
#define A01    v22
#define A11    v23
#define A02    v24
#define A12    v25
#define A03    v26
#define A13    v27

#define C00    v20
#define C01    v21
#define C02    v22
#define C03    v23
#define C04    v24
#define C05    v25
#define C06    v26
#define C07    v27
#define C08    v28
#define C09    v29

// aliases of C00, ..., C09 for easier readability
#define C10    v20
#define C11    v21
#define C12    v22
#define C13    v23
#define C14    v24
#define C15    v25
#define C16    v26
#define C17    v27
#define C18    v28
#define C19    v29

#define AB00   v0
#define AB10   v1
#define AB01   v2
#define AB11   v3
#define AB02   v4
#define AB12   v5
#define AB03   v6
#define AB13   v7
#define AB04   v8
#define AB14   v9
#define AB05   v10
#define AB15   v11
#define AB06   v12
#define AB16   v13
#define AB07   v14
#define AB17   v15
#define AB08   v16
#define AB18   v17
#define AB09   v18
#define AB19   v19

#define rs_c   a6
#define cs_c   a7

REALNAME:
    #include "rviv_save_registers.h"

    vsetvli s0, zero, VTYPE, m1, ta, ma
    csrr s0, vlenb

    // Zero-initialize accumulators
    vxor.vv AB00, AB00, AB00
    vxor.vv AB10, AB10, AB10
    vxor.vv AB01, AB01, AB01
    vxor.vv AB11, AB11, AB11
    vxor.vv AB02, AB02, AB02
    vxor.vv AB12, AB12, AB12
    vxor.vv AB03, AB03, AB03
    vxor.vv AB13, AB13, AB13
    vxor.vv AB04, AB04, AB04
    vxor.vv AB14, AB14, AB14
    vxor.vv AB05, AB05, AB05
    vxor.vv AB15, AB15, AB15
    vxor.vv AB06, AB06, AB06
    vxor.vv AB16, AB16, AB16
    vxor.vv AB07, AB07, AB07
    vxor.vv AB17, AB17, AB17
    vxor.vv AB07, AB07, AB07
    vxor.vv AB17, AB17, AB17
    vxor.vv AB08, AB08, AB08
    vxor.vv AB18, AB18, AB18
    vxor.vv AB09, AB09, AB09
    vxor.vv AB19, AB19, AB19

    // Handle k == 0
    beqz loop_counter, MULTIPLYBETA

    // Set up pointers to A
    add A10_ptr, A00_ptr, s0
    slli s0, s0, 1 // length of a column of A in bytes
    add A01_ptr, A00_ptr, s0
    add A11_ptr, A10_ptr, s0

    li tmp, 3
    ble loop_counter, tmp, TAIL_UNROLL_2

    // Set up the remaining pointer to A
    add A02_ptr, A01_ptr, s0
    add A12_ptr, A11_ptr, s0
    add A03_ptr, A02_ptr, s0
    add A13_ptr, A12_ptr, s0

    // Preload A and B
    // Load A(:,l)
    VLE A00, (A00_ptr)
    VLE A10, (A10_ptr)

    // Load B(l,0:9)
    FLOAD B00, 0*DATASIZE(B_row_ptr)
    FLOAD B01, 1*DATASIZE(B_row_ptr)
    FLOAD B02, 2*DATASIZE(B_row_ptr)
    FLOAD B03, 3*DATASIZE(B_row_ptr)
    FLOAD B04, 4*DATASIZE(B_row_ptr)
    FLOAD B05, 5*DATASIZE(B_row_ptr)
    FLOAD B06, 6*DATASIZE(B_row_ptr)
    FLOAD B07, 7*DATASIZE(B_row_ptr)
    FLOAD B08, 8*DATASIZE(B_row_ptr)
    FLOAD B09, 9*DATASIZE(B_row_ptr)

LOOP_UNROLL_4:
    addi loop_counter, loop_counter, -4

    vfmacc.vf AB00, B00, A00   // AB(:,0:3) += A(:,l) * B(l,0:3)
    vfmacc.vf AB01, B01, A00
    vfmacc.vf AB02, B02, A00
    vfmacc.vf AB03, B03, A00
    vfmacc.vf AB10, B00, A10
    vfmacc.vf AB11, B01, A10
    vfmacc.vf AB12, B02, A10
    vfmacc.vf AB13, B03, A10

    // Load B(l+1,0:9)
    FLOAD B10, 10*DATASIZE(B_row_ptr)
    FLOAD B11, 11*DATASIZE(B_row_ptr)
    FLOAD B12, 12*DATASIZE(B_row_ptr)
    FLOAD B13, 13*DATASIZE(B_row_ptr)
    FLOAD B14, 14*DATASIZE(B_row_ptr)
    FLOAD B15, 15*DATASIZE(B_row_ptr)
    FLOAD B16, 16*DATASIZE(B_row_ptr)
    FLOAD B17, 17*DATASIZE(B_row_ptr)
    FLOAD B18, 18*DATASIZE(B_row_ptr)
    FLOAD B19, 19*DATASIZE(B_row_ptr)
    addi B_row_ptr, B_row_ptr, 20*DATASIZE

    // Load A(:,l+1)
    VLE A01, (A01_ptr)
    VLE A11, (A11_ptr)

    vfmacc.vf AB04, B04, A00   // AB(:,4:9) += A(:,l) * B(l,4:9)
    vfmacc.vf AB05, B05, A00
    vfmacc.vf AB06, B06, A00
    vfmacc.vf AB07, B07, A00
    vfmacc.vf AB08, B08, A00
    vfmacc.vf AB09, B09, A00
    vfmacc.vf AB14, B04, A10
    vfmacc.vf AB15, B05, A10
    vfmacc.vf AB16, B06, A10
    vfmacc.vf AB17, B07, A10
    vfmacc.vf AB18, B08, A10
    vfmacc.vf AB19, B09, A10

    // Load A(:,l+2)
    VLE A02, (A02_ptr)
    VLE A12, (A12_ptr)

    vfmacc.vf AB00, B10, A01   // AB(:,0:3) += A(:,l+1) * B(l+1,0:3)
    vfmacc.vf AB01, B11, A01
    vfmacc.vf AB02, B12, A01
    vfmacc.vf AB03, B13, A01
    vfmacc.vf AB10, B10, A11
    vfmacc.vf AB11, B11, A11
    vfmacc.vf AB12, B12, A11
    vfmacc.vf AB13, B13, A11

    // Load B(l+2,0:9)
    FLOAD B20, 0*DATASIZE(B_row_ptr)
    FLOAD B21, 1*DATASIZE(B_row_ptr)
    FLOAD B22, 2*DATASIZE(B_row_ptr)
    FLOAD B23, 3*DATASIZE(B_row_ptr)
    FLOAD B24, 4*DATASIZE(B_row_ptr)
    FLOAD B25, 5*DATASIZE(B_row_ptr)
    FLOAD B26, 6*DATASIZE(B_row_ptr)
    FLOAD B27, 7*DATASIZE(B_row_ptr)
    FLOAD B28, 8*DATASIZE(B_row_ptr)
    FLOAD B29, 9*DATASIZE(B_row_ptr)

    vfmacc.vf AB04, B14, A01   // AB(:,4:9) += A(:,l+1) * B(l+1,4:9)
    vfmacc.vf AB05, B15, A01
    vfmacc.vf AB06, B16, A01
    vfmacc.vf AB07, B17, A01
    vfmacc.vf AB08, B18, A01
    vfmacc.vf AB09, B19, A01
    vfmacc.vf AB14, B14, A11
    vfmacc.vf AB15, B15, A11
    vfmacc.vf AB16, B16, A11
    vfmacc.vf AB17, B17, A11
    vfmacc.vf AB18, B18, A11
    vfmacc.vf AB19, B19, A11

    // Load B(l+3,0:9)
    FLOAD B30, 10*DATASIZE(B_row_ptr)
    FLOAD B31, 11*DATASIZE(B_row_ptr)
    FLOAD B32, 12*DATASIZE(B_row_ptr)
    FLOAD B33, 13*DATASIZE(B_row_ptr)
    FLOAD B34, 14*DATASIZE(B_row_ptr)
    FLOAD B35, 15*DATASIZE(B_row_ptr)
    FLOAD B36, 16*DATASIZE(B_row_ptr)
    FLOAD B37, 17*DATASIZE(B_row_ptr)
    FLOAD B38, 18*DATASIZE(B_row_ptr)
    FLOAD B39, 19*DATASIZE(B_row_ptr)
    addi B_row_ptr, B_row_ptr, 20*DATASIZE

    vfmacc.vf AB00, B20, A02   // AB(:,0:3) += A(:,l+2) * B(l+2,0:3)
    vfmacc.vf AB01, B21, A02
    vfmacc.vf AB02, B22, A02
    vfmacc.vf AB03, B23, A02
    vfmacc.vf AB10, B20, A12
    vfmacc.vf AB11, B21, A12
    vfmacc.vf AB12, B22, A12
    vfmacc.vf AB13, B23, A12

    // Load A(:,l+3)
    VLE A03, (A03_ptr)
    VLE A13, (A13_ptr)

    vfmacc.vf AB04, B24, A02   // AB(:,4:9) += A(:,l+2) * B(l+2,4:9)
    vfmacc.vf AB05, B25, A02
    vfmacc.vf AB06, B26, A02
    vfmacc.vf AB07, B27, A02
    vfmacc.vf AB08, B28, A02
    vfmacc.vf AB09, B29, A02
    vfmacc.vf AB14, B24, A12
    vfmacc.vf AB15, B25, A12
    vfmacc.vf AB16, B26, A12
    vfmacc.vf AB17, B27, A12
    vfmacc.vf AB18, B28, A12
    vfmacc.vf AB19, B29, A12

    // Advance pointers to A(:,l+5:l+6)
    add A00_ptr, A03_ptr, s0
    add A10_ptr, A13_ptr, s0
    add A01_ptr, A00_ptr, s0
    add A11_ptr, A10_ptr, s0

    vfmacc.vf AB00, B30, A03   // AB(:,0:3) += A(:,l+3) * B(l+3,0:3)
    vfmacc.vf AB01, B31, A03
    vfmacc.vf AB02, B32, A03
    vfmacc.vf AB03, B33, A03
    vfmacc.vf AB10, B30, A13
    vfmacc.vf AB11, B31, A13
    vfmacc.vf AB12, B32, A13
    vfmacc.vf AB13, B33, A13

    vfmacc.vf AB04, B34, A03   // AB(:,4:9) += A(:,l+3) * B(l+3,4:9)
    vfmacc.vf AB05, B35, A03
    vfmacc.vf AB06, B36, A03
    vfmacc.vf AB07, B37, A03
    vfmacc.vf AB08, B38, A03
    vfmacc.vf AB09, B39, A03
    vfmacc.vf AB14, B34, A13
    vfmacc.vf AB15, B35, A13
    vfmacc.vf AB16, B36, A13
    vfmacc.vf AB17, B37, A13
    vfmacc.vf AB18, B38, A13
    vfmacc.vf AB19, B39, A13

    li tmp, 3
    ble loop_counter, tmp, TAIL_UNROLL_2

    // Load A and B for the next iteration
    // Load A(:,l+5)
    VLE A00, (A00_ptr)
    VLE A10, (A10_ptr)

    // Load B(l+5,0:9)
    FLOAD B00, 0*DATASIZE(B_row_ptr)
    FLOAD B01, 1*DATASIZE(B_row_ptr)
    FLOAD B02, 2*DATASIZE(B_row_ptr)
    FLOAD B03, 3*DATASIZE(B_row_ptr)
    FLOAD B04, 4*DATASIZE(B_row_ptr)
    FLOAD B05, 5*DATASIZE(B_row_ptr)
    FLOAD B06, 6*DATASIZE(B_row_ptr)
    FLOAD B07, 7*DATASIZE(B_row_ptr)
    FLOAD B08, 8*DATASIZE(B_row_ptr)
    FLOAD B09, 9*DATASIZE(B_row_ptr)

    // Advance pointers to A(:,l+7:l+8)
    add A02_ptr, A01_ptr, s0
    add A12_ptr, A11_ptr, s0
    add A03_ptr, A02_ptr, s0
    add A13_ptr, A12_ptr, s0

    j LOOP_UNROLL_4

TAIL_UNROLL_2: // loop_counter <= 3
    li tmp, 1
    ble loop_counter, tmp, TAIL_UNROLL_1

    addi loop_counter, loop_counter, -2

    // Load B(l,0:9)
    FLOAD B00, 0*DATASIZE(B_row_ptr)
    FLOAD B01, 1*DATASIZE(B_row_ptr)
    FLOAD B02, 2*DATASIZE(B_row_ptr)
    FLOAD B03, 3*DATASIZE(B_row_ptr)
    FLOAD B04, 4*DATASIZE(B_row_ptr)
    FLOAD B05, 5*DATASIZE(B_row_ptr)
    FLOAD B06, 6*DATASIZE(B_row_ptr)
    FLOAD B07, 7*DATASIZE(B_row_ptr)
    FLOAD B08, 8*DATASIZE(B_row_ptr)
    FLOAD B09, 9*DATASIZE(B_row_ptr)

    // Load A(:,l)
    VLE A00, (A00_ptr)
    VLE A10, (A10_ptr)

    vfmacc.vf AB00, B00, A00   // AB(:,0:3) += A(:,l) * B(l,0:3)
    vfmacc.vf AB01, B01, A00
    vfmacc.vf AB02, B02, A00
    vfmacc.vf AB03, B03, A00
    vfmacc.vf AB10, B00, A10
    vfmacc.vf AB11, B01, A10
    vfmacc.vf AB12, B02, A10
    vfmacc.vf AB13, B03, A10

    // Load B(l+1,0:9)
    FLOAD B10, 10*DATASIZE(B_row_ptr)
    FLOAD B11, 11*DATASIZE(B_row_ptr)
    FLOAD B12, 12*DATASIZE(B_row_ptr)
    FLOAD B13, 13*DATASIZE(B_row_ptr)
    FLOAD B14, 14*DATASIZE(B_row_ptr)
    FLOAD B15, 15*DATASIZE(B_row_ptr)
    FLOAD B16, 16*DATASIZE(B_row_ptr)
    FLOAD B17, 17*DATASIZE(B_row_ptr)
    FLOAD B18, 18*DATASIZE(B_row_ptr)
    FLOAD B19, 19*DATASIZE(B_row_ptr)
    addi B_row_ptr, B_row_ptr, 20*DATASIZE

    vfmacc.vf AB04, B04, A00   // AB(:,4:9) += A(:,l) * B(l,4:9)
    vfmacc.vf AB05, B05, A00
    vfmacc.vf AB06, B06, A00
    vfmacc.vf AB07, B07, A00
    vfmacc.vf AB08, B08, A00
    vfmacc.vf AB09, B09, A00
    vfmacc.vf AB14, B04, A10
    vfmacc.vf AB15, B05, A10
    vfmacc.vf AB16, B06, A10
    vfmacc.vf AB17, B07, A10
    vfmacc.vf AB18, B08, A10
    vfmacc.vf AB19, B09, A10

    // Load A(:,l+1)
    VLE A01, (A01_ptr)
    VLE A11, (A11_ptr)

    vfmacc.vf AB00, B10, A01   // AB(:,0:3) += A(:,l+1) * B(l+1,0:3)
    vfmacc.vf AB01, B11, A01
    vfmacc.vf AB02, B12, A01
    vfmacc.vf AB03, B13, A01
    vfmacc.vf AB10, B10, A11
    vfmacc.vf AB11, B11, A11
    vfmacc.vf AB12, B12, A11
    vfmacc.vf AB13, B13, A11

    // Advance pointers to A(:,l+2)
    add A00_ptr, A01_ptr, s0
    add A10_ptr, A11_ptr, s0

    vfmacc.vf AB04, B14, A01   // AB(:,4:9) += A(:,l+1) * B(l+1,4:9)
    vfmacc.vf AB05, B15, A01
    vfmacc.vf AB06, B16, A01
    vfmacc.vf AB07, B17, A01
    vfmacc.vf AB08, B18, A01
    vfmacc.vf AB09, B19, A01
    vfmacc.vf AB14, B14, A11
    vfmacc.vf AB15, B15, A11
    vfmacc.vf AB16, B16, A11
    vfmacc.vf AB17, B17, A11
    vfmacc.vf AB18, B18, A11
    vfmacc.vf AB19, B19, A11

    li tmp, 1
    ble loop_counter, tmp, TAIL_UNROLL_1

TAIL_UNROLL_1: // loop_counter <= 1
    beqz loop_counter, MULTIPLYALPHA

    // Load B(l,0:9)
    FLOAD B00, 0*DATASIZE(B_row_ptr)
    FLOAD B01, 1*DATASIZE(B_row_ptr)
    FLOAD B02, 2*DATASIZE(B_row_ptr)
    FLOAD B03, 3*DATASIZE(B_row_ptr)
    FLOAD B04, 4*DATASIZE(B_row_ptr)
    FLOAD B05, 5*DATASIZE(B_row_ptr)
    FLOAD B06, 6*DATASIZE(B_row_ptr)
    FLOAD B07, 7*DATASIZE(B_row_ptr)
    FLOAD B08, 8*DATASIZE(B_row_ptr)
    FLOAD B09, 9*DATASIZE(B_row_ptr)

    // Load A(:,l)
    VLE A00, (A00_ptr)
    VLE A10, (A10_ptr)

    vfmacc.vf AB00, B00, A00   // AB(:,0:3) += A(:,l) * B(l,0:3)
    vfmacc.vf AB01, B01, A00
    vfmacc.vf AB02, B02, A00
    vfmacc.vf AB03, B03, A00
    vfmacc.vf AB10, B00, A10
    vfmacc.vf AB11, B01, A10
    vfmacc.vf AB12, B02, A10
    vfmacc.vf AB13, B03, A10

    vfmacc.vf AB04, B04, A00   // AB(:,4:9) += A(:,l) * B(l,4:9)
    vfmacc.vf AB05, B05, A00
    vfmacc.vf AB06, B06, A00
    vfmacc.vf AB07, B07, A00
    vfmacc.vf AB08, B08, A00
    vfmacc.vf AB09, B09, A00
    vfmacc.vf AB14, B04, A10
    vfmacc.vf AB15, B05, A10
    vfmacc.vf AB16, B06, A10
    vfmacc.vf AB17, B07, A10
    vfmacc.vf AB18, B08, A10
    vfmacc.vf AB19, B09, A10

MULTIPLYALPHA:
    FLOAD ALPHA, (a1)

    // Multiply with alpha
    vfmul.vf AB00, AB00, ALPHA
    vfmul.vf AB01, AB01, ALPHA
    vfmul.vf AB02, AB02, ALPHA
    vfmul.vf AB03, AB03, ALPHA
    vfmul.vf AB04, AB04, ALPHA
    vfmul.vf AB05, AB05, ALPHA
    vfmul.vf AB06, AB06, ALPHA
    vfmul.vf AB07, AB07, ALPHA
    vfmul.vf AB08, AB08, ALPHA
    vfmul.vf AB09, AB09, ALPHA

    vfmul.vf AB10, AB10, ALPHA
    vfmul.vf AB11, AB11, ALPHA
    vfmul.vf AB12, AB12, ALPHA
    vfmul.vf AB13, AB13, ALPHA
    vfmul.vf AB14, AB14, ALPHA
    vfmul.vf AB15, AB15, ALPHA
    vfmul.vf AB16, AB16, ALPHA
    vfmul.vf AB17, AB17, ALPHA
    vfmul.vf AB18, AB18, ALPHA
    vfmul.vf AB19, AB19, ALPHA

MULTIPLYBETA:
    FLOAD BETA,  (a4)
    FZERO(fzero)

    // Set up pointers to columns of C
    add C1_ptr, C0_ptr, cs_c
    add C2_ptr, C1_ptr, cs_c
    add C3_ptr, C2_ptr, cs_c
    add C4_ptr, C3_ptr, cs_c
    add C5_ptr, C4_ptr, cs_c
    add C6_ptr, C5_ptr, cs_c
    add C7_ptr, C6_ptr, cs_c
    add C8_ptr, C7_ptr, cs_c
    add C9_ptr, C8_ptr, cs_c

    FEQ tmp, BETA, fzero
    beq tmp, zero, BETANOTZERO

BETAZERO:
    VSE AB00, (C0_ptr)
    VSE AB01, (C1_ptr)
    VSE AB02, (C2_ptr)
    VSE AB03, (C3_ptr)
    VSE AB04, (C4_ptr)
    VSE AB05, (C5_ptr)
    VSE AB06, (C6_ptr)
    VSE AB07, (C7_ptr)
    VSE AB08, (C8_ptr)
    VSE AB09, (C9_ptr)

    // Advance columns pointer by VLEN
    add C0_ptr, C0_ptr, rs_c
    add C1_ptr, C1_ptr, rs_c
    add C2_ptr, C2_ptr, rs_c
    add C3_ptr, C3_ptr, rs_c
    add C4_ptr, C4_ptr, rs_c
    add C5_ptr, C5_ptr, rs_c
    add C6_ptr, C6_ptr, rs_c
    add C7_ptr, C7_ptr, rs_c
    add C8_ptr, C8_ptr, rs_c
    add C9_ptr, C9_ptr, rs_c

    VSE AB10, (C0_ptr)
    VSE AB11, (C1_ptr)
    VSE AB12, (C2_ptr)
    VSE AB13, (C3_ptr)
    VSE AB14, (C4_ptr)
    VSE AB15, (C5_ptr)
    VSE AB16, (C6_ptr)
    VSE AB17, (C7_ptr)
    VSE AB18, (C8_ptr)
    VSE AB19, (C9_ptr)

    j END

BETANOTZERO:
    VLE C00, (C0_ptr)
    VLE C01, (C1_ptr)
    VLE C02, (C2_ptr)
    VLE C03, (C3_ptr)
    VLE C04, (C4_ptr)
    VLE C05, (C5_ptr)
    VLE C06, (C6_ptr)
    VLE C07, (C7_ptr)
    VLE C08, (C8_ptr)
    VLE C09, (C9_ptr)

    vfmacc.vf AB00, BETA, C00
    vfmacc.vf AB01, BETA, C01
    vfmacc.vf AB02, BETA, C02
    vfmacc.vf AB03, BETA, C03

    VSE AB00, (C0_ptr)
    VSE AB01, (C1_ptr)
    VSE AB02, (C2_ptr)
    VSE AB03, (C3_ptr)

    add C0_ptr, C0_ptr, rs_c   // Advance pointers by VLEN
    add C1_ptr, C1_ptr, rs_c
    add C2_ptr, C2_ptr, rs_c
    add C3_ptr, C3_ptr, rs_c

    vfmacc.vf AB04, BETA, C04
    vfmacc.vf AB05, BETA, C05
    vfmacc.vf AB06, BETA, C06
    vfmacc.vf AB07, BETA, C07
    vfmacc.vf AB08, BETA, C08
    vfmacc.vf AB09, BETA, C09

    VSE AB04, (C4_ptr)
    VSE AB05, (C5_ptr)
    VSE AB06, (C6_ptr)
    VSE AB07, (C7_ptr)
    VSE AB08, (C8_ptr)
    VSE AB09, (C9_ptr)

    add C4_ptr, C4_ptr, rs_c
    add C5_ptr, C5_ptr, rs_c
    add C6_ptr, C6_ptr, rs_c
    add C7_ptr, C7_ptr, rs_c
    add C8_ptr, C8_ptr, rs_c
    add C9_ptr, C9_ptr, rs_c

    VLE C10, (C0_ptr)
    VLE C11, (C1_ptr)
    VLE C12, (C2_ptr)
    VLE C13, (C3_ptr)
    VLE C14, (C4_ptr)
    VLE C15, (C5_ptr)
    VLE C16, (C6_ptr)
    VLE C17, (C7_ptr)
    VLE C18, (C8_ptr)
    VLE C19, (C9_ptr)

    vfmacc.vf AB10, BETA, C10
    vfmacc.vf AB11, BETA, C11
    vfmacc.vf AB12, BETA, C12
    vfmacc.vf AB13, BETA, C13
    vfmacc.vf AB14, BETA, C14
    vfmacc.vf AB15, BETA, C15
    vfmacc.vf AB16, BETA, C16
    vfmacc.vf AB17, BETA, C17
    vfmacc.vf AB18, BETA, C18
    vfmacc.vf AB19, BETA, C19

    VSE AB10, (C0_ptr)
    VSE AB11, (C1_ptr)
    VSE AB12, (C2_ptr)
    VSE AB13, (C3_ptr)
    VSE AB14, (C4_ptr)
    VSE AB15, (C5_ptr)
    VSE AB16, (C6_ptr)
    VSE AB17, (C7_ptr)
    VSE AB18, (C8_ptr)
    VSE AB19, (C9_ptr)

END:
    #include "rviv_restore_registers.h"
    ret

