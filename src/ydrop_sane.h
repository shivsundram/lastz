//-------+---------+---------+---------+---------+---------+---------+--------=
//
// File: ydrop_sane.h
//
//----------
//
// Two textbook ports of lastz's ydrop_one_sided_align(). Both are
// bit-identical to it on the v0 scope (forward extension only, no
// leftSeg/rightSeg, no active-segment masking, trimToPeak == true).
//
//   1. ydrop_one_sided_align_impl_sane
//      "Wasteful but obvious" reference: three full (M+1)x(N+1) score
//      matrices and one full (M+1)x(N+1) link byte tape, all per-call.
//      ~13 * (M+1)*(N+1) bytes. Easiest to read alongside the algebraic
//      recurrences.
//
//   2. ydrop_one_sided_align_impl_sane_double_buffered
//      Same algorithm, lastz-matching memory profile: two (N+1)-cell
//      sweep-row buffers each for C and D, scalar I, and a band-compact
//      link byte tape indexed via a per-row offset table. Per-call
//      memory is O(N) for scores and O(sum of band widths) for links.
//
// Both are tested against lastz in bench/test_ydrop.c.
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

// Reference impl: full (M+1)x(N+1) score matrices and link tape.
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

// Lastz-matching memory: sweep-row C/D buffers + band-compact link tape.
score ydrop_one_sided_align_impl_sane_double_buffered (
    const u8   *A,         unspos M,
    const u8   *B,         unspos N,
    scorerow   *allSub,
    score       gapOpen,
    score       gapExtend,
    score       yDrop,
    editscript **script_out,
    unspos     *end1_out,
    unspos     *end2_out);

//----------
// Per-thread reusable scratch buffers for the pooled variant of the
// double-buffered impl. Avoids the ~8 malloc/free pairs per call (× 2
// halves per anchor × thousands of anchors) of the unpooled version.
//
// Lifecycle:
//   - ydrop_scratch_new()    : allocate an empty scratch (no buffers yet).
//   - ydrop_scratch_free()   : release all backing memory.
//   - ydrop_one_sided_align_impl_sane_double_buffered_pooled() takes a
//     scratch by reference and grows the buffers in-place as needed.
//     Subsequent calls reuse what's already there; first call sizes
//     everything once.
//
// Thread safety:
//   Each ydrop_scratch is single-owner. Callers driving parallel y-drop
//   must use one scratch per thread (e.g. an array indexed by
//   omp_get_thread_num).
//----------

typedef struct ydrop_scratch ydrop_scratch;

ydrop_scratch *ydrop_scratch_new  (void);
void           ydrop_scratch_free (ydrop_scratch *s);

// Identical observable behavior to ydrop_one_sided_align_impl_sane_double_buffered
// but reuses the buffers in *scratch across calls. scratch must be non-NULL.
score ydrop_one_sided_align_impl_sane_double_buffered_pooled (
    const u8       *A,         unspos M,
    const u8       *B,         unspos N,
    scorerow       *allSub,
    score           gapOpen,
    score           gapExtend,
    score           yDrop,
    editscript    **script_out,
    unspos         *end1_out,
    unspos         *end2_out,
    ydrop_scratch  *scratch);

#ifdef __cplusplus
}
#endif

#endif // ydrop_sane_H
