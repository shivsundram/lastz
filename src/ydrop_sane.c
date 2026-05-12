//-------+---------+---------+---------+---------+---------+---------+--------=
//
// File: ydrop_sane.c
//
//----------
//
// Two textbook ports of lastz's ydrop_one_sided_align(), both bit-identical
// to it on the v0 scope (forward extension only, no leftSeg/rightSeg, no
// active-segment masking, trimToPeak == true).
//
//   1. ydrop_one_sided_align_impl_sane
//      The "wasteful but obvious" reference: three full (M+1)x(N+1)
//      matrices for C, D, I plus one full (M+1)x(N+1) link byte tape.
//      Memory ~13 * (M+1)*(N+1) bytes. The clearest possible mapping
//      from the algebraic recurrences to code; intended as the reading
//      copy when you want to understand what the GPU kernel must compute.
//
//   2. ydrop_one_sided_align_impl_sane_double_buffered
//      The lastz-matching memory profile: two (N+1)-cell sweep-row
//      buffers each for C and D, scalar I, and a band-compact link byte
//      tape indexed via a per-row offset table (cell (r,c)'s link byte
//      lives at lnkTape[lnkRow[r] + c]). Per-call memory is O(N) for
//      scores and O(sum of band widths) for links. The algorithm is
//      identical to (1) — only the storage layout differs.
//
// Both functions are validated bit-identical against
// ydrop_one_sided_align in bench/test_ydrop.c.
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
//     (i.e., when the cell's link is cFromC).
//   - In the "we can improve C" branch, both iExtend and dExtend are
//     forced ON unconditionally.
//   - In the "we cannot improve C" branch, dExtend is set iff
//     (D[r][c] - gapE) >= (C[r][c] - gapOE), and similarly for iExtend.
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

//==========================================================================
//
// ydrop_one_sided_align_impl_sane
//
//   Reference / "textbook" implementation: three full (M+1)x(N+1) score
//   matrices for C, D, I plus one full (M+1)x(N+1) link byte tape. No
//   sweep-row aliasing, no dp/dq pointer trick, no static state, no
//   module globals. Per-call malloc_or_die and free_if_valid.
//
//   Memory at score_type=I (4 B/cell scores, 1 B/cell link):
//     13 * (M+1) * (N+1) bytes  =  ~13 MB at M=N=1000, ~52 MB at M=N=2000.
//   The companion test driver caps inputs at <= 2000 bp.
//
//   For the lastz-matching memory profile (sweep-row scores +
//   band-compact link tape), use the _double_buffered variant below.
//
//   Same convention as gapped_extend.c::ydrop_one_sided_align: A[0] and
//   B[0] are the anchor base; DP runs on A[1..M] vs B[1..N]. Cell (r,c)
//   aligns A[1..r] with B[1..c]. Returns the best score reached;
//   *end1_out/*end2_out give that score's location (1-based offsets
//   past the anchor). The edit script is APPENDED to *script_out.
//
//==========================================================================

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


//==========================================================================
//
// ydrop_one_sided_align_impl_sane_double_buffered
//
//   Same algorithm as ydrop_one_sided_align_impl_sane above, but with
//   lastz-matching memory layout (see gapped_extend.c note 3 and line
//   3775):
//
//   - Sweep-row score buffers: two (N+1)-cell arrays each for C and D,
//     swapped between rows. I is a scalar carried through the inner loop.
//   - Band-compact link tape: one byte tape, indexed via a per-row
//     offset table. Cell (r, c)'s link byte is at lnkTape[lnkRow[r] + c].
//     Only cells in the live band [LY[r], RY[r]) consume tape bytes;
//     cells outside the band do not exist in memory.
//   - LY/RY arrays of length M+1.
//
//   Per-call memory:
//     - Score buffers:    4 * (N+1) * sizeof(score)   (~16 KB at N=1000).
//     - Link tape:        O(sum_{r} band_width). For typical y-drop runs
//                         (band ~500 cells), ~500 KB at M=1000, ~5 MB at
//                         M=10000.
//     - lnkRow/LY/RY:     O(M).
//
//   Unlike lastz, we do per-call malloc/free (no module-static reuse),
//   keeping the impl state-free at the cost of some allocator churn.
//
//   The band-compact tape indexing trick (mirrors gapped_extend.c:3775):
//
//     lnkRow[r] = (tape offset where row r's segment begins) - LY[r]
//
//     Then for c in [LY[r], RY[r]):
//       lnkTape[lnkRow[r] + c]
//     lands at byte (c - LY[r]) of row r's segment.
//
//     When LY[r] > current_offset, the stored lnkRow[r] is a very large
//     unsigned, but adding c (which is >= LY[r]) brings it back into the
//     valid range modulo 2^64. The traceback only ever indexes with c in
//     [LY[r], RY[r]), so we never address out of bounds.
//
//==========================================================================

// Initial size for the link tape and growth factor on overflow.
#define LNKTAPE_INITIAL  (64u * 1024u)
#define LNKTAPE_GROWTH   2u

score ydrop_one_sided_align_impl_sane_double_buffered (
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

    // yDropTail = how many insertions in a row before -yDrop fires.
    // Used as per-row tape headroom so realloc rarely needs to fire
    // during a row.
    score yDropTail = (gapE > 0 ? yDrop / gapE : 0) + 16;

    // ---- Sweep-row score buffers (4 * (N+1) cells total) ----
    size_t row_cells = (size_t)N + 1;
    score *C_prev = (score*) malloc_or_die ("ydrop_sane_db C_prev",
                                            row_cells * sizeof(score));
    score *C_curr = (score*) malloc_or_die ("ydrop_sane_db C_curr",
                                            row_cells * sizeof(score));
    score *D_prev = (score*) malloc_or_die ("ydrop_sane_db D_prev",
                                            row_cells * sizeof(score));
    score *D_curr = (score*) malloc_or_die ("ydrop_sane_db D_curr",
                                            row_cells * sizeof(score));

    for (size_t k = 0; k < row_cells; k++)
        { C_prev[k] = C_curr[k] = D_prev[k] = D_curr[k] = negInf; }

    // ---- Band-compact link tape ----
    size_t  lnkTapeSize = LNKTAPE_INITIAL;
    u8     *lnkTape     = (u8*) malloc_or_die ("ydrop_sane_db lnkTape",
                                               lnkTapeSize);
    size_t *lnkRow      = (size_t*) malloc_or_die ("ydrop_sane_db lnkRow",
                                                   (size_t)(M + 2) * sizeof(size_t));
    size_t  lnkPos      = 0;

    // ---- Per-row band bounds ----
    unspos *LY = (unspos*) malloc_or_die ("ydrop_sane_db LY",
                                          (size_t)(M + 1) * sizeof(unspos));
    unspos *RY = (unspos*) malloc_or_die ("ydrop_sane_db RY",
                                          (size_t)(M + 1) * sizeof(unspos));

    score  bestScore = 0;
    unspos end1 = 0, end2 = 0;

    // ---- Row 0: pure-insertion prefix until y-drop bites ----
    LY[0]     = 0;
    lnkRow[0] = lnkPos - (size_t)LY[0];

    if (lnkPos + (size_t)N + 1 > lnkTapeSize)
        {
        size_t newSize = lnkTapeSize;
        while (newSize < lnkPos + (size_t)N + 1) newSize *= LNKTAPE_GROWTH;
        lnkTape = (u8*) realloc_or_die ("ydrop_sane_db lnkTape",
                                        lnkTape, newSize);
        lnkTapeSize = newSize;
        }

    C_prev[0]              = 0;
    lnkTape[lnkRow[0] + 0] = 0;
    lnkPos++;

    unspos rowZeroStop = 1;
    {
    score prev_c = 0;
    score next_c = -gapOE;
    for (unspos col = 1; col <= N && prev_c >= -yDrop; col++)
        {
        C_prev[col]              = next_c;
        lnkTape[lnkRow[0] + col] = cFromI;
        lnkPos++;
        rowZeroStop = col + 1;
        prev_c      = next_c;
        next_c     -= gapE;
        }
    }
    RY[0] = rowZeroStop;

    // ---- Main sweep: row-major ----
    unspos r;
    for (r = 1; r <= M; r++)
        {
        unspos initialLY = LY[r-1];
        unspos rowRight  = RY[r-1];

        // Ensure tape has room for this row's worst-case writes.
        size_t tbNeeded = (size_t)(rowRight - initialLY)
                        + (size_t)yDropTail
                        + 1;
        if (lnkPos + tbNeeded > lnkTapeSize)
            {
            size_t newSize = lnkTapeSize;
            while (newSize < lnkPos + tbNeeded) newSize *= LNKTAPE_GROWTH;
            lnkTape = (u8*) realloc_or_die ("ydrop_sane_db lnkTape",
                                            lnkTape, newSize);
            lnkTapeSize = newSize;
            }

        LY[r]     = initialLY;
        lnkRow[r] = lnkPos - (size_t)initialLY;

        score *sub_row  = allSub[A[r]];
        long   npCol    = (long)LY[r] - 1;
        score  I_scalar = negInf;

        for (unspos c = LY[r]; c < rowRight; c++)
            {
            score Drc = MAX (C_prev[c] - gapOE,
                             D_prev[c] - gapE);

            score Irc = (c == 0)
                      ? negInf
                      : MAX (C_curr[c-1] - gapOE,
                             I_scalar  - gapE);

            score diag = (c == 0)
                       ? negInf
                       : C_prev[c-1] + sub_row[B[c]];

            score Crc;  u8 link;
            int   pruned = 0;

            if (Drc > diag || Irc > diag)
                {
                if (Drc >= Irc) { Crc = Drc; link = cFromD | iExtend | dExtend; }
                else            { Crc = Irc; link = cFromI | iExtend | dExtend; }
                if (Crc < bestScore - yDrop) pruned = 1;
                }
            else
                {
                Crc = diag;  link = cFromC;

                if (Crc < bestScore - yDrop)
                    pruned = 1;
                else
                    {
                    if (Crc >= bestScore)
                        { bestScore = Crc; end1 = r; end2 = c; }
                    score cOpen = Crc - gapOE;
                    if ((Drc - gapE) >= cOpen) link |= dExtend;
                    if ((Irc - gapE) >= cOpen) link |= iExtend;
                    }
                }

            if (pruned)
                {
                C_curr[c] = negInf;
                D_curr[c] = negInf;
                I_scalar  = negInf;
                lnkTape[lnkRow[r] + c] = 0;
                if (c == LY[r]) LY[r]++;
                }
            else
                {
                C_curr[c] = Crc;
                D_curr[c] = Drc;
                I_scalar  = Irc;
                lnkTape[lnkRow[r] + c] = link;
                npCol     = (long)c;
                }

            lnkPos++;
            }

        if ((long)rowRight > npCol + 1)
            {
            RY[r] = (unspos)(npCol + 1);
            }
        else
            {
            RY[r] = rowRight;
            unspos col = rowRight;
            if (col >= 1 && col <= N)
                {
                score i_next = MAX (C_curr[col-1] - gapOE,
                                    I_scalar     - gapE);
                while ((i_next >= bestScore - yDrop) && (col <= N))
                    {
                    C_curr[col]              = i_next;
                    I_scalar                 = i_next;
                    D_curr[col]              = i_next - gapOE;
                    lnkTape[lnkRow[r] + col] = cFromI;
                    lnkPos++;
                    RY[r]                    = col + 1;
                    col++;
                    if (col > N) break;
                    i_next -= gapE;
                    }
                }
            }

        // Right-boundary termination: grow the band by one extra column
        // if there's room. The boundary cell stays negInf in C_curr/D_curr;
        // no tape write — traceback never reads this cell.
        if (RY[r] <= N)
            {
            RY[r]++;
            C_curr[RY[r] - 1] = negInf;
            D_curr[RY[r] - 1] = negInf;
            }

        // Ghost the cell immediately left of the live band, so that the
        // next row's diag read at c = LY[r+1] (which equals LY[r] before
        // any shrink) sees negInf instead of a stale value left over from
        // an older row.
        if (LY[r] >= 1)
            {
            C_curr[LY[r] - 1] = negInf;
            D_curr[LY[r] - 1] = negInf;
            }

        if (LY[r] >= RY[r]) break;

        // Swap sweep-row buffers for the next row.
        {
        score *tmp;
        tmp = C_prev; C_prev = C_curr; C_curr = tmp;
        tmp = D_prev; D_prev = D_curr; D_curr = tmp;
        }
        }

    // ---- Traceback: walk lnkTape backward from (end1, end2) ----
    {
    unspos row = end1, col = end2;
    u8 prevOp = 0, op;
    while (row >= 1 || col > 0)
        {
        u8 link = lnkTape[lnkRow[row] + col];
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

    free_if_valid ("ydrop_sane_db RY",      RY);
    free_if_valid ("ydrop_sane_db LY",      LY);
    free_if_valid ("ydrop_sane_db lnkRow",  lnkRow);
    free_if_valid ("ydrop_sane_db lnkTape", lnkTape);
    free_if_valid ("ydrop_sane_db D_curr",  D_curr);
    free_if_valid ("ydrop_sane_db D_prev",  D_prev);
    free_if_valid ("ydrop_sane_db C_curr",  C_curr);
    free_if_valid ("ydrop_sane_db C_prev",  C_prev);

    return bestScore;
    }
