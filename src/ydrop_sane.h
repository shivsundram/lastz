//-------+---------+---------+---------+---------+---------+---------+--------=
//
// File: ydrop_sane.h
//
//----------
//
// ydrop_sane--
//   Textbook 2-D row-major Gotoh y-drop, intended as a readable reference
//   implementation of ydrop_one_sided_align(). Uses three full (M+1)x(N+1)
//   matrices for C, D, I and one (M+1)x(N+1) link byte tape, all allocated
//   per-call via malloc_or_die and freed before return.
//
//   v0 scope: forward extension only, no leftSeg/rightSeg, no active-segment
//   masking, trimToPeak == true. The intent is bit-identical agreement with
//   lastz's ydrop_one_sided_align on calls that fall in this scope (which
//   includes all the test cases in bench/test_ydrop.c).
//
//----------

#ifndef ydrop_sane_H
#define ydrop_sane_H

#include "build_options.h"
#include "utilities.h"
#include "dna_utilities.h"
#include "edit_script.h"

#ifdef __cplusplus
extern "C" {
#endif

score ydrop_one_sided_align_impl_sane (
    const u8   *A,         unspos M,
    const u8   *B,         unspos N,
    scorerow   *allSub,
    score       gapOpen,
    score       gapExtend,
    score       yDrop,
    editscript **script_out,
    unspos     *end1_out,
    unspos     *end2_out);

#ifdef __cplusplus
}
#endif

#endif // ydrop_sane_H
