//-------+---------+---------+---------+---------+---------+---------+--------=
//
// File: ydrop_sane.c
//
//----------
//
// Textbook 2-D row-major Gotoh y-drop, intended as a readable reference
// implementation of ydrop_one_sided_align() in gapped_extend.c.
//
// Three full (M+1)x(N+1) matrices for C, D, I, and one (M+1)x(N+1) link
// byte tape. No sweep-row aliasing, no dp/dq pointer trick, no static
// state, no module globals. Per-call malloc_or_die and free_if_valid
// (the same allocators lastz uses), so this participates in lastz's leak
// tracking when that's enabled.
//
// Memory at score_type=I (4 B/cell scores, 1 B/cell link):
//   13 * (M+1) * (N+1) bytes  =  ~13 MB at M=N=1000, ~52 MB at M=N=2000.
// The companion test driver caps inputs at <= 2000 bp.
//
// Recurrences (standard Gotoh affine-gap, see note 5 in gapped_extend.c):
//   D[r][c] = max( C[r-1][c]   - gapOE,  D[r-1][c]   - gapE )
//   I[r][c] = max( C[r][c-1]   - gapOE,  I[r][c-1]   - gapE )
//   C[r][c] = max( C[r-1][c-1] + sub[A[r]][B[c]],  D[r][c],  I[r][c] )
//
// Y-drop pruning: a cell with C[r][c] < bestScore - yDrop is set to negInf
// (with link cleared) and the per-row [LY[r], RY[r]) feasible range is
// shrunk; subsequent rows inherit this band.
//
// Link byte format (matches gapped_extend.c note 7):
//   bits 0-1: cFromC (0) | cFromI (1) | cFromD (2)
//   bit 2:    iExtend
//   bit 3:    dExtend
// The bestScore endpoint and the iExtend/dExtend semantics match lastz's
// ydrop_one_sided_align exactly:
//   - bestScore is updated ONLY in the "we cannot improve C" branch
//     (i.e., when the cell's link is cFromC). Cells whose final C value
//     came from D or I are NOT eligible to set the alignment endpoint.
//   - In the "we can improve C" branch, both iExtend and dExtend are
//     forced ON unconditionally (matches lastz's branch comment).
//   - In the "we cannot improve C" branch, dExtend is set iff
//     (D[r][c] - gapE) >= (C[r][c] - gapOE), i.e., extending the
//     existing D into the cell BELOW beats opening a new D from the
//     current C. iExtend is set iff (I[r][c] - gapE) >= (C[r][c] - gapOE)
//     for the cell to the RIGHT.
//
//----------

#include <stdlib.h>
#include <string.h>

#include "build_options.h"
#include "utilities.h"
#include "dna_utilities.h"
#include "edit_script.h"

#include "ydrop_sane.h"

// Traceback link bits — same encoding as gapped_extend.c.
#define cFromC   0
#define cFromI   1
#define cFromD   2
#define iExtend  4
#define dExtend  8
#define cidBits  (cFromC | cFromI | cFromD)

// Match gapped_extend.c's local alias.
#define negInf negInfinity

#ifndef MAX
#define MAX(a,b) (((a) >= (b)) ? (a) : (b))
#endif

//----------
//
// ydrop_one_sided_align_impl_sane--
//
// Forward (rightward + downward) gapped extension from the anchor at the
// origin of (A, B). Matches lastz's ydrop_one_sided_align convention:
//   - A[0] and B[0] are the ANCHOR base, NOT part of the DP. (Caller
//     would typically pass &seq1[anchor1], &seq2[anchor2], plus the
//     remaining lengths after the anchor.)
//   - The DP processes A[1..M] vs B[1..N]. Cell [r][c] is the alignment
//     of A[1..r] with B[1..c], so the substitution score for that cell
//     is sub[A[r]][B[c]].
//
// Returns the best score reached. *end1_out / *end2_out give the position
// of that score (1-based offsets into the post-anchor extension; so
// end1=k means we extended k bases on A past the anchor). The edit
// script is APPENDED to *script_out.
//
//----------

score ydrop_one_sided_align_impl_sane (
    const u8   *A,         unspos M,
    const u8   *B,         unspos N,
    scorerow   *allSub,
    score       gapOpen,
    score       gapExtend,
    score       yDrop,
    editscript **script_out,
    unspos     *end1_out,
    unspos     *end2_out)
    {
    score gapE  = gapExtend;
    score gapOE = gapOpen + gapExtend;

    if (M == 0 || N == 0)
        { *end1_out = *end2_out = 0; return 0; }

    // ---- Allocate the four full (M+1)x(N+1) matrices (per-call) ----
    size_t cells = (size_t)(M + 1) * (size_t)(N + 1);
    score *C    = (score*) malloc_or_die ("ydrop_sane C",    cells * sizeof(score));
    score *D    = (score*) malloc_or_die ("ydrop_sane D",    cells * sizeof(score));
    score *I    = (score*) malloc_or_die ("ydrop_sane I",    cells * sizeof(score));
    u8    *lnk  = (u8*)    malloc_or_die ("ydrop_sane link", cells * sizeof(u8));

    // Convenience: row-major flat-array index. NB: depends on N being captured.
    #define IDX(r,c) ((size_t)(r) * (size_t)(N + 1) + (size_t)(c))

    // Initialize everything to "unreachable" sentinel.
    for (size_t k = 0; k < cells; k++)
        { C[k] = D[k] = I[k] = negInf; lnk[k] = 0; }

    // ---- Row 0: pure-insertion prefix until y-drop bites ----
    // Lastz convention (gapped_extend.c lines 3677-3697):
    //   for (col=1; col<=N && cTemp >= -yDrop; col++) { set C[0][col] = c; ... }
    // i.e., we keep filling as long as the PREVIOUS cell was still above
    // -yDrop. The fill therefore extends one cell PAST the cutoff.
    score  bestScore = 0;
    unspos end1 = 0, end2 = 0;

    C[IDX(0,0)] = 0;
    lnk[IDX(0,0)] = 0;

    unspos rowZeroStop = 1;   // one-past-rightmost initialized col on row 0
    {
    score prev_c = 0;          // C[0][0]
    score next_c = -gapOE;     // candidate C[0][1]
    for (unspos col = 1; col <= N && prev_c >= -yDrop; col++)
        {
        C[IDX(0,col)]   = next_c;
        lnk[IDX(0,col)] = cFromI;
        rowZeroStop = col + 1;
        prev_c = next_c;
        next_c -= gapE;
        }
    }

    // Per-row feasible bounds.
    //   LY[r] = leftmost live col on row r
    //   RY[r] = one past rightmost live col on row r
    unspos *LY = (unspos*) malloc_or_die ("ydrop_sane LY", (M + 1) * sizeof(unspos));
    unspos *RY = (unspos*) malloc_or_die ("ydrop_sane RY", (M + 1) * sizeof(unspos));
    LY[0] = 0;
    RY[0] = rowZeroStop;

    // ---- Main DP sweep: row-major ----
    // Convention for each row r:
    //   - LY[r] starts at LY[r-1] and may shrink right via consecutive
    //     left-edge prunes (lastz: `if (col == LY) LY++;`).
    //   - The inner loop iterates c from LY[r] to RY[r-1]-1. Mid-row prunes
    //     write negInf into C/D/I but do NOT break the loop or shrink the
    //     band immediately. We track npCol = highest c that was NOT pruned.
    //   - After the inner loop:
    //       if (RY[r-1] > npCol+1)  →  y-drop fired before reaching the right
    //                                 edge. Shrink: RY[r] = npCol+1.
    //       else                    →  reached the right edge. Prolong this
    //                                 row by a chain of pure insertions.
    unspos r;
    for (r = 1; r <= M; r++)
        {
        LY[r] = LY[r-1];
        unspos rowRight = RY[r-1];

        // Substitution-score row vector for target base A[r]. Lastz
        // convention: A[0] is the anchor, A[1..M] is the extension; cell
        // [r][c] aligns A[r] with B[c].
        score *sub_row = allSub[A[r]];

        // npCol = last column on this row where the cell was NOT pruned.
        // Use a signed sentinel so we can express "no live cell yet".
        long npCol = (long) LY[r] - 1;

        for (unspos c = LY[r]; c < rowRight; c++)
            {
            // (1) Affine recurrences.
            score Drc = MAX (C[IDX(r-1,c)] - gapOE,
                             D[IDX(r-1,c)] - gapE);

            score Irc = (c == 0)
                      ? negInf
                      : MAX (C[IDX(r,c-1)] - gapOE,
                             I[IDX(r,c-1)] - gapE);

            score diag = (c == 0)
                       ? negInf
                       : C[IDX(r-1,c-1)] + sub_row[B[c]];

            // (2) C: 3-way max + traceback link, matching lastz's branching.
            score Crc;  u8 link;
            int   pruned = 0;

            if (Drc > diag || Irc > diag)
                {
                // "We CAN improve C" branch.
                if (Drc >= Irc) { Crc = Drc; link = cFromD | iExtend | dExtend; }
                else            { Crc = Irc; link = cFromI | iExtend | dExtend; }

                // y-drop prune (no bestScore update on this branch).
                if (Crc < bestScore - yDrop) pruned = 1;
                }
            else
                {
                // "We CANNOT improve C" branch — C stays as the diagonal value.
                Crc = diag;  link = cFromC;

                if (Crc < bestScore - yDrop)
                    pruned = 1;
                else
                    {
                    // bestScore update — ONLY on the cFromC branch (matches lastz).
                    if (Crc >= bestScore)
                        { bestScore = Crc; end1 = r; end2 = c; }

                    // Conditional iExtend / dExtend bits.
                    score cOpen = Crc - gapOE;
                    if ((Drc - gapE) >= cOpen) link |= dExtend;
                    if ((Irc - gapE) >= cOpen) link |= iExtend;
                    }
                }

            if (pruned)
                {
                C[IDX(r,c)]   = negInf;
                D[IDX(r,c)]   = negInf;
                I[IDX(r,c)]   = negInf;
                lnk[IDX(r,c)] = 0;
                if (c == LY[r]) LY[r]++;     // chain-shrink the left edge
                // (no break: keep iterating; the next cell may still be alive)
                }
            else
                {
                C[IDX(r,c)]   = Crc;
                D[IDX(r,c)]   = Drc;
                I[IDX(r,c)]   = Irc;
                lnk[IDX(r,c)] = link;
                npCol         = (long) c;
                }
            }

        if ((long) rowRight > npCol + 1)
            {
            // y-drop fired before reaching the previous row's right edge.
            RY[r] = (unspos) (npCol + 1);
            }
        else
            {
            // Reached the right edge without y-drop firing. Prolong this
            // row to the right via a chain of pure insertions while
            //   I[r][col] = max(C[r][col-1] - gapOE, I[r][col-1] - gapE)
            // remains >= bestScore - yDrop and col <= N.
            // (lastz lines 3919-3937)
            RY[r]      = rowRight;
            unspos col = rowRight;
            if (col >= 1 && col <= N)
                {
                score i_next = MAX (C[IDX(r,col-1)] - gapOE,
                                    I[IDX(r,col-1)] - gapE);
                while ((i_next >= bestScore - yDrop) && (col <= N))
                    {
                    C[IDX(r,col)]   = i_next;
                    I[IDX(r,col)]   = i_next;
                    D[IDX(r,col)]   = i_next - gapOE;
                    lnk[IDX(r,col)] = cFromI;
                    RY[r]           = col + 1;
                    col++;
                    if (col > N) break;
                    i_next -= gapE;
                    }
                }
            }

        // Right-boundary termination (lastz lines 3943-3951): grow the
        // band by one extra column if there's room, with that cell set to
        // negInf. This makes the next row's inner loop reach one column
        // past the previous live cell, where it may find a fresh diagonal
        // match. The cell itself stays negInf (already initialized so).
        if (RY[r] <= N) RY[r]++;

        if (LY[r] >= RY[r]) break;
        }

    // ---- Traceback: walk lnk[] backward from (end1, end2) ----
    {
    unspos row = end1, col = end2;
    u8 prevOp = 0, op;
    while (row >= 1 || col > 0)
        {
        u8 link = lnk[IDX(row,col)];
        op = link & cidBits;
        if ((prevOp == cFromI) && ((link & iExtend) != 0)) op = cFromI;
        if ((prevOp == cFromD) && ((link & dExtend) != 0)) op = cFromD;

        if      (op == cFromI) {           col--; edit_script_ins (script_out, 1); }
        else if (op == cFromD) { row--;           edit_script_del (script_out, 1); }
        else                   { row--;    col--; edit_script_sub (script_out, 1); }
        prevOp = op;
        }
    }

    *end1_out = end1;
    *end2_out = end2;

    #undef IDX

    free_if_valid ("ydrop_sane RY",   RY);
    free_if_valid ("ydrop_sane LY",   LY);
    free_if_valid ("ydrop_sane link", lnk);
    free_if_valid ("ydrop_sane I",    I);
    free_if_valid ("ydrop_sane D",    D);
    free_if_valid ("ydrop_sane C",    C);

    return bestScore;
    }
