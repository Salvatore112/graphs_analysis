//------------------------------------------------------------------------------
// LAGraph_msf.c
//------------------------------------------------------------------------------

// LAGraph, (c) 2019-2022 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause
//
// For additional details (including references to third party source code and
// other files) see the LICENSE file or contact permission@sei.cmu.edu. See
// Contributors.txt for a full list of contributors. Created, in part, with
// funding and support from the U.S. Government (see Acknowledgments.txt file).
// DM22-0790

// Contributed by Yongzhe Zhang (zyz915@gmail.com)

//------------------------------------------------------------------------------

/**
 * Code is based on Boruvka's minimum spanning forest algorithm
 */

// FIXME: is this ready for src?  It uses global values, so not yet ready.
// FIXME: Reduce_assign is slow.  See src/algorithm/LG_CC_FastSV6/7.

#ifndef LG_INTERNAL_H
#define LG_INTERNAL_H

//------------------------------------------------------------------------------
// include files
//------------------------------------------------------------------------------

#include "LAGraph.h"
#include <ctype.h>
#undef I

#if defined(__linux__)
#include <malloc.h>
#endif

#define LG_RESTRICT LAGRAPH_RESTRICT

//------------------------------------------------------------------------------
// string macros
//------------------------------------------------------------------------------

#define LG_XSTR(x) LG_STR(x)
#define LG_STR(x) #x

//------------------------------------------------------------------------------
// string matching
//------------------------------------------------------------------------------

#define MATCH(s1, s2, n) (strncmp(s1, s2, n) == 0)
#define MATCHNAME(s1, s2) MATCH(s1, s2, LAGRAPH_MAX_NAME_LEN)

//------------------------------------------------------------------------------
// typedefs
//------------------------------------------------------------------------------

// LG_void: used in place of (void *), but valid for pointer arithmetic
typedef unsigned char LG_void;

//------------------------------------------------------------------------------
// LG_CLEAR_MSG: clear the error msg string
//------------------------------------------------------------------------------

// When an LAGraph method starts, it first clears the caller's msg string.
#define LG_CLEAR_MSG                                                                                                   \
    {                                                                                                                  \
        if (msg != NULL)                                                                                               \
            msg[0] = '\0';                                                                                             \
    }

//------------------------------------------------------------------------------
// LG_ERROR_MSG: set the error msg string
//------------------------------------------------------------------------------

// When an LAGraph method encounters an error, it can report details in the
// msg.  This is normally done via LG_ASSERT_MSG.  For example:

/*
    if (src < 0 || src >= n)
    {
        LG_ERROR_MSG ("Source node %ld must be in range 0 to n-1, "
            "where n = %ld is the number of nodes in the graph.", src, n) ;
        return (GrB_INVALID_INDEX) ;
    }
    // or, with a simpler message:
    LG_ASSERT_MSG (src >= 0 && src < n, GrB_INVALID_INDEX, "invalid src node") ;
*/

#define LG_ERROR_MSG(...)                                                                                              \
    {                                                                                                                  \
        if (msg != NULL && msg[0] == '\0') {                                                                           \
            snprintf(msg, LAGRAPH_MSG_LEN, __VA_ARGS__);                                                               \
        }                                                                                                              \
    }

//------------------------------------------------------------------------------
// LG_FREE_WORK: free all workspace
//------------------------------------------------------------------------------

#ifndef LG_FREE_WORK
#define LG_FREE_WORK ;
#endif

//------------------------------------------------------------------------------
// LG_FREE_ALL: free all workspace and all output arguments, on error
//------------------------------------------------------------------------------

#ifndef LG_FREE_ALL
#define LG_FREE_ALL                                                                                                    \
    {                                                                                                                  \
        LG_FREE_WORK;                                                                                                  \
    }
#endif

//------------------------------------------------------------------------------
// GRB_CATCH: catch an error from GraphBLAS
//------------------------------------------------------------------------------

// A simple GRB_CATCH macro to be used by GRB_TRY.  If an LAGraph function
// wants something else, then #define a GRB_CATCH macro before the #include
// "LG_internal.h" statement.

#ifndef GRB_CATCH
#define GRB_CATCH(info)                                                                                                \
    {                                                                                                                  \
        LG_ERROR_MSG("GraphBLAS failure (file %s, line %d): info: %d", __FILE__, __LINE__, info);                      \
        LG_FREE_ALL;                                                                                                   \
        return (info);                                                                                                 \
    }
#endif

//------------------------------------------------------------------------------
// LAGRAPH_CATCH: catch an error from LAGraph
//------------------------------------------------------------------------------

// A simple LAGRAPH_CATCH macro to be used by LAGRAPH_TRY.  If an LAGraph
// function wants something else, then #define a LAGRAPH_CATCH macro before the
// #include "LG_internal.h" statement.

#ifndef LAGRAPH_CATCH
#define LAGRAPH_CATCH(status)                                                                                          \
    {                                                                                                                  \
        LG_ERROR_MSG("LAGraph failure (file %s, line %d): status: %d", __FILE__, __LINE__, status);                    \
        LG_FREE_ALL;                                                                                                   \
        return (status);                                                                                               \
    }
#endif

//------------------------------------------------------------------------------
// LG_ASSERT_MSGF: assert an expression is true, and return if it is false
//------------------------------------------------------------------------------

// Identical to LG_ASSERT_MSG, except this allows a printf-style formatted
// message.

#define LG_ASSERT_MSGF(expression, error_status, expression_format, ...)                                               \
    {                                                                                                                  \
        if (!(expression)) {                                                                                           \
            LG_ERROR_MSG("LAGraph failure (file %s, line %d): " expression_format, __FILE__, __LINE__, __VA_ARGS__);   \
            LG_FREE_ALL;                                                                                               \
            return (error_status);                                                                                     \
        }                                                                                                              \
    }

//------------------------------------------------------------------------------
// LG_ASSERT_MSG: assert an expression is true, and return if it is false
//------------------------------------------------------------------------------

// Identical to LG_ASSERT, except this allows a different string to be
// included in the message.

#define LG_ASSERT_MSG(expression, error_status, expression_message)                                                    \
    LG_ASSERT_MSGF(expression, error_status, "%s", expression_message)

//------------------------------------------------------------------------------
// LG_ASSERT: assert an expression is true, and return if it is false
//------------------------------------------------------------------------------

// LAGraph methods can use this assertion macro for simple errors.

#define LG_ASSERT(expression, error_status)                                                                            \
    {                                                                                                                  \
        if (!(expression)) {                                                                                           \
            LG_ERROR_MSG("LAGraph assertion \"%s\" failed (file %s, line %d):"                                         \
                         " status: %d",                                                                                \
                         LG_XSTR(expression), __FILE__, __LINE__, error_status);                                       \
            LG_FREE_ALL;                                                                                               \
            return (error_status);                                                                                     \
        }                                                                                                              \
    }

//------------------------------------------------------------------------------
// LG_TRY: check a condition and return on error
//------------------------------------------------------------------------------

// The msg is not modified.  This should be used when an LAGraph method calls
// another one.

#define LG_TRY(LAGraph_method)                                                                                         \
    {                                                                                                                  \
        int LAGraph_status = LAGraph_method;                                                                           \
        if (LAGraph_status < 0) {                                                                                      \
            LG_FREE_ALL;                                                                                               \
            return (LAGraph_status);                                                                                   \
        }                                                                                                              \
    }

//------------------------------------------------------------------------------
// LG_CLEAR_MSG_AND_BASIC_ASSERT: clear msg and do basic tests of a graph
//------------------------------------------------------------------------------

#define LG_CLEAR_MSG_AND_BASIC_ASSERT(G, msg)                                                                          \
    {                                                                                                                  \
        LG_CLEAR_MSG;                                                                                                  \
        LG_ASSERT(G != NULL, GrB_NULL_POINTER);                                                                        \
        LG_ASSERT_MSG(G->A != NULL, LAGRAPH_INVALID_GRAPH, "graph adjacency matrix is NULL");                          \
        LG_ASSERT_MSG(G->kind >= LAGraph_ADJACENCY_UNDIRECTED && G->kind <= LAGraph_ADJACENCY_DIRECTED,                \
                      LAGRAPH_INVALID_GRAPH, "graph kind invalid");                                                    \
    }

//------------------------------------------------------------------------------
// FPRINTF: fprintf and check result
//------------------------------------------------------------------------------

#define FPRINTF(f, ...)                                                                                                \
    {                                                                                                                  \
        LG_ASSERT_MSG(fprintf(f, __VA_ARGS__) >= 0, LAGRAPH_IO_ERROR, "Unable to write to file");                      \
    }

//------------------------------------------------------------------------------
// code development settings
//------------------------------------------------------------------------------

// turn off debugging; do not edit these three lines
#ifndef NDEBUG
#define NDEBUG
#endif

// These flags are used for code development.  Uncomment them as needed.

// to turn on debugging, uncomment this line:
// #undef NDEBUG

#undef ASSERT

#ifndef NDEBUG

// debugging enabled
#ifdef MATLAB_MEX_FILE
// debugging when LAGraph is part of a mexFunction
#define ASSERT(x)                                                                                                      \
    {                                                                                                                  \
        if (!(x))                                                                                                      \
            mexErrMsgTxt("failure: " __FILE__ " line: ");                                                              \
    }
#else
#include <assert.h>
#define ASSERT(x) assert(x);
#endif

#else

// debugging disabled
#define ASSERT(x)

#endif

//------------------------------------------------------------------------------
// LG_Multiply_size_t:  c = a*b but check for overflow
//------------------------------------------------------------------------------

static bool LG_Multiply_size_t // true if ok, false if overflow
    (size_t *c,                // c = a*b, or zero if overflow occurs
     const size_t a, const size_t b) {

    ASSERT(c != NULL);

    (*c) = 0;
    if (a == 0 || b == 0) {
        return (true);
    }

    if (a > SIZE_MAX / 2 || b > SIZE_MAX / 2) {
        // a or b are out of range
        return (false);
    }

    // a + b is now safe to compute
    if ((a + b) > (SIZE_MAX / LAGRAPH_MIN(a, b))) {
        // a * b may overflow
        return (false);
    }

    // a * b will not overflow
    (*c) = a * b;
    return (true);
}

//------------------------------------------------------------------------------
// Matrix Market format
//------------------------------------------------------------------------------

// %%MatrixMarket matrix <fmt> <type> <storage> uses the following enums:

typedef enum {
    MM_coordinate,
    MM_array,
} MM_fmt_enum;

typedef enum { MM_real, MM_integer, MM_complex, MM_pattern } MM_type_enum;

typedef enum { MM_general, MM_symmetric, MM_skew_symmetric, MM_hermitian } MM_storage_enum;

// maximum length of each line in the Matrix Market file format

// The MatrixMarket format specificies a maximum line length of 1024.
// This is currently sufficient for GraphBLAS but will need to be relaxed
// if this function is extended to handle arbitrary user-defined types.
#define MMLEN 1024
#define MAXLINE MMLEN + 6

//------------------------------------------------------------------------------
// LG_PART and LG_PARTITION: definitions for partitioning an index range
//------------------------------------------------------------------------------

LAGRAPH_PUBLIC extern int LG_nthreads_outer; // # of threads to use at the higher level of a nested
                                             // parallel region in LAGraph.  Default: 1.

LAGRAPH_PUBLIC extern int LG_nthreads_inner; // # of threads to use at the lower level of a nested
                                             // parallel region, or to use inside GraphBLAS.
                                             // Default: the value obtained by omp_get_max_threads
                                             // if OpenMP is in use, or 1 otherwise.

// LG_PART and LG_PARTITION:  divide the index range 0:n-1 uniformly
// for nthreads.  LG_PART(tid,n,nthreads) is the first index for thread tid.
#define LG_PART(tid, n, nthreads) (((tid) * ((double)(n))) / ((double)(nthreads)))

// thread tid will operate on the range k1:(k2-1)
#define LG_PARTITION(k1, k2, n, tid, nthreads)                                                                         \
    k1 = ((tid) == 0) ? 0 : LG_PART((tid), n, nthreads);                                                               \
    k2 = ((tid) == (nthreads) - 1) ? (n) : LG_PART((tid) + 1, n, nthreads)

//------------------------------------------------------------------------------
// LG_eslice: uniform partition of e items to each task
//------------------------------------------------------------------------------

static inline void LG_eslice(int64_t *Slice,  // array of size ntasks+1
                             int64_t e,       // number items to partition amongst the tasks
                             const int ntasks // # of tasks
) {
    Slice[0] = 0;
    for (int tid = 0; tid < ntasks; tid++) {
        Slice[tid] = LG_PART(tid, e, ntasks);
    }
    Slice[ntasks] = e;
}

//------------------------------------------------------------------------------
// definitions for sorting functions
//------------------------------------------------------------------------------

// All of the LG_qsort_* functions are single-threaded, by design.  The
// LG_msort* functions are parallel.  None of these sorting methods are
// guaranteed to be stable.  These functions are contributed by Tim Davis, and
// are derived from SuiteSparse:GraphBLAS.  Functions named LG_* are not
// meant to be accessible by end users of LAGraph.

#define LG_BASECASE (64 * 1024)

//------------------------------------------------------------------------------
// LG_msort1: sort array of size n
//------------------------------------------------------------------------------

// LG_msort1 sorts an int64_t array of size n in ascending order.

int LG_msort1(
    // input/output:
    int64_t *A_0, // size n array
    // input:
    const int64_t n, char *msg);

//------------------------------------------------------------------------------
// LG_msort2: sort two arrays of size n
//------------------------------------------------------------------------------

// LG_msort2 sorts two int64_t arrays A of size n in ascending order.
// The arrays are kept in the same order, where the pair (A_0 [k], A_1 [k]) is
// treated as a single pair.  The pairs are sorted by the first value A_0,
// with ties broken by A_1.

int LG_msort2(
    // input/output:
    int64_t *A_0, // size n array
    int64_t *A_1, // size n array
    // input:
    const int64_t n, char *msg);

//------------------------------------------------------------------------------
// LG_msort3: sort three arrays of size n
//------------------------------------------------------------------------------

// LG_msort3 sorts three int64_t arrays A of size n in ascending order.
// The arrays are kept in the same order, where the triplet (A_0 [k], A_1 [k],
// A_2 [k]) is treated as a single triplet.  The triplets are sorted by the
// first value A_0, with ties broken by A_1, and then by A_2 if the values of
// A_0 and A_1 are identical.

int LG_msort3(
    // input/output:
    int64_t *A_0, // size n array
    int64_t *A_1, // size n array
    int64_t *A_2, // size n array
    // input:
    const int64_t n, char *msg);

void LG_qsort_1a               // sort array A of size 1-by-n
    (int64_t *LG_RESTRICT A_0, // size n array
     const int64_t n);

void LG_qsort_2                // sort array A of size 2-by-n, using 2 keys (A [0:1][])
    (int64_t *LG_RESTRICT A_0, // size n array
     int64_t *LG_RESTRICT A_1, // size n array
     const int64_t n);

void LG_qsort_3                // sort array A of size 3-by-n, using 3 keys (A [0:2][])
    (int64_t *LG_RESTRICT A_0, // size n array
     int64_t *LG_RESTRICT A_1, // size n array
     int64_t *LG_RESTRICT A_2, // size n array
     const int64_t n);

//------------------------------------------------------------------------------
// LG_lt_1: sorting comparator function, one key
//------------------------------------------------------------------------------

// A [a] and B [b] are keys of one integer.

// LG_lt_1 returns true if A [a] < B [b], for LG_qsort_1b

#define LG_lt_1(A_0, a, B_0, b) (A_0[a] < B_0[b])

//------------------------------------------------------------------------------
// LG_lt_2: sorting comparator function, two keys
//------------------------------------------------------------------------------

// A [a] and B [b] are keys of two integers.

// LG_lt_2 returns true if A [a] < B [b], for LG_qsort_2 and LG_msort_2b

#define LG_lt_2(A_0, A_1, a, B_0, B_1, b)                                                                              \
    ((A_0[a] < B_0[b]) ? (true)                                                                                        \
                       : ((A_0[a] == B_0[b]) ? (/* primary key is the same; tie-break on the 2nd key */                \
                                                (A_1[a] < B_1[b]))                                                     \
                                             : (false)))

//------------------------------------------------------------------------------
// LG_lt_3: sorting comparator function, three keys
//------------------------------------------------------------------------------

// A [a] and B [b] are keys of three integers.

// LG_lt_3 returns true if A [a] < B [b], for LG_qsort_3 and LG_msort_3b

#define LG_lt_3(A_0, A_1, A_2, a, B_0, B_1, B_2, b)                                                                    \
    ((A_0[a] < B_0[b]) ? (true)                                                                                        \
                       : ((A_0[a] == B_0[b]) ? (/* primary key is the same; tie-break on the 2nd and 3rd key */        \
                                                LG_lt_2(A_1, A_2, a, B_1, B_2, b))                                     \
                                             : (false)))

//------------------------------------------------------------------------------
// LG_eq_*: sorting comparator function, three keys
//------------------------------------------------------------------------------

// A [a] and B [b] are keys of two or three integers.
// LG_eq_* returns true if A [a] == B [b]

#define LG_eq_3(A_0, A_1, A_2, a, B_0, B_1, B_2, b) ((A_0[a] == B_0[b]) && (A_1[a] == B_1[b]) && (A_2[a] == B_2[b]))

#define LG_eq_2(A_0, A_1, a, B_0, B_1, b) ((A_0[a] == B_0[b]) && (A_1[a] == B_1[b]))

#define LG_eq_1(A_0, a, B_0, b) ((A_0[a] == B_0[b]))

//------------------------------------------------------------------------------
// count entries on the diagonal of a matrix
//------------------------------------------------------------------------------

int LG_nself_edges(
    // output
    int64_t *nself_edges, // # of entries
    // input
    GrB_Matrix A, // matrix to count
    char *msg     // error message
);

//------------------------------------------------------------------------------
// simple and portable random number generator (internal use only)
//------------------------------------------------------------------------------

#define LG_RANDOM15_MAX 32767
#define LG_RANDOM60_MAX ((1ULL << 60) - 1)

// return a random number between 0 and LG_RANDOM15_MAX
GrB_Index LG_Random15(uint64_t *seed);

// return a random uint64_t, in range 0 to LG_RANDOM60_MAX
GrB_Index LG_Random60(uint64_t *seed);

//------------------------------------------------------------------------------
// LG_KindName: return the name of a kind
//------------------------------------------------------------------------------

// LG_KindName: return the name of a graph kind.  For example, if given
// LAGraph_ADJACENCY_UNDIRECTED, the string "undirected" is returned.

int LG_KindName(
    // output:
    char *name, // name of the kind (user provided array of size at least
                // LAGRAPH_MAX_NAME_LEN)
    // input:
    LAGraph_Kind kind, // graph kind
    char *msg);

//------------------------------------------------------------------------------

// # of entries to print for LAGraph_Matrix_Print and LAGraph_Vector_Print
#define LG_SHORT_LEN 30

//------------------------------------------------------------------------------
// GrB_get/set for SuiteSparse extensions
//------------------------------------------------------------------------------

#if LAGRAPH_SUITESPARSE

// SuiteSparse:GraphBLAS v8.1 or later
#define LG_HYPERSPARSE GxB_HYPERSPARSE
#define LG_SPARSE GxB_SPARSE
#define LG_BITMAP GxB_BITMAP
#define LG_FULL GxB_FULL
#define LG_SET_FORMAT_HINT(object, sparsity) GrB_set(object, (int32_t)(sparsity), GxB_SPARSITY_CONTROL)
#define LG_GET_FORMAT_HINT(A, status) GrB_get(A, (int32_t *)status, GxB_SPARSITY_STATUS)
#define LG_SET_NTHREADS(nthreads) GrB_set(GrB_GLOBAL, (int32_t)(nthreads), GxB_NTHREADS)
#define LG_SET_HYPER_SWITCH(A, hyper) GrB_set(A, hyper, GxB_HYPER_SWITCH)
#define LG_GET_HYPER_SWITCH(A, hyper) GrB_get(A, hyper, GxB_HYPER_SWITCH)
#define LG_SET_BURBLE(burble) GrB_set(GrB_GLOBAL, (int32_t)(burble), GxB_BURBLE)
#define LG_GET_LIBRARY_DATE(date) GrB_get(GrB_GLOBAL, (char *)date, GxB_LIBRARY_DATE)

#if defined(GRAPHBLAS_HAS_CUDA)
// the LG_brutal_malloc family of methods.
#define LG_BRUTAL_TESTS 0
#else
#define LG_BRUTAL_TESTS 1
#endif

#else

// vanilla GraphBLAS
#define LG_HYPERSPARSE 1
#define LG_SPARSE 2
#define LG_BITMAP 4
#define LG_FULL 8
#define LG_SET_FORMAT_HINT(object, sparsity) GrB_SUCCESS
#define LG_SET_NTHREADS(nthreads) GrB_SUCCESS
#define LG_SET_HYPER_SWITCH(A, hyper) GrB_SUCCESS
#define LG_GET_HYPER_SWITCH(A, hyper) GrB_SUCCESS
#define LG_GET_FORMAT_HINT(A, status) GrB_SUCCESS
#define LG_SET_BURBLE(burble) GrB_SUCCESS
#define LG_GET_LIBRARY_DATE(date) GrB_SUCCESS
#define LG_BRUTAL_TESTS 0

#endif

#endif

#include <LAGraph.h>
#include <LAGraphX.h>

//****************************************************************************
// encode each edge into a single uint64_t
static void combine(void *z, const void *x, const void *y) {
    *(uint64_t *)z = ((*(uint64_t *)x) << 32) + (*(uint64_t *)y);
}

static void get_fst(void *y, const void *x) { *(uint64_t *)y = (*(uint64_t *)x) >> 32; }

static void get_snd(void *y, const void *x) { *(uint64_t *)y = (*(uint64_t *)x) & INT_MAX; }

//****************************************************************************
// FIXME: Reduce_assign is slow.  See src/algorithm/LG_CC_FastSV6.

#undef LG_FREE_ALL
#define LG_FREE_ALL LAGraph_Free((void **)&mem, msg);

// w[index[i]] = min(w[index[i]], s[i]) for i in [0..n-1]
static GrB_Info Reduce_assign(GrB_Vector w, GrB_Vector s, GrB_Index *index, GrB_Index n, char *msg) {
    GrB_Index *mem = NULL;
    LG_TRY(LAGraph_Malloc((void **)&mem, n * 3, sizeof(GrB_Index), msg));
    GrB_Index *ind = mem, *sval = mem + n, *wval = sval + n;
    LG_TRY(GrB_Vector_extractTuples(ind, wval, &n, w));
    LG_TRY(GrB_Vector_extractTuples(ind, sval, &n, s));
    for (GrB_Index i = 0; i < n; i++)
        if (sval[i] < wval[index[i]])
            wval[index[i]] = sval[i];
    LG_TRY(GrB_Vector_clear(w));
    LG_TRY(GrB_Vector_build(w, ind, wval, n, GrB_PLUS_UINT64));
    LG_FREE_ALL;
    return GrB_SUCCESS;
}

//****************************************************************************
// global C arrays (for implementing various GrB_IndexUnaryOps)
static GrB_Index *weight = NULL, *parent = NULL, *partner = NULL;

// generate solution:
// for each element A(i, j), it is selected if
//   1. weight[i] == A(i, j)    -- where weight[i] stores i's minimum edge weight
//   2. parent[j] == partner[i] -- j belongs to the specified connected component

void f1(bool *z, const void *x, GrB_Index i, GrB_Index j, const void *thunk) {
    uint64_t *aij = (uint64_t *)x;
    (*z) = (weight[i] == *aij) && (parent[j] == partner[i]);
}

// edge removal:
// A(i, j) is removed when parent[i] == parent[j]

void f2(bool *z, const void *x, GrB_Index i, GrB_Index j, const void *thunk) { (*z) = (parent[i] != parent[j]); }

//****************************************************************************

#undef LG_FREE_ALL
#define LG_FREE_ALL                                                                                                    \
    {                                                                                                                  \
        GrB_free(&S);                                                                                                  \
        GrB_free(&T);                                                                                                  \
        LAGraph_Free((void **)&I, msg);                                                                                \
        LAGraph_Free((void **)&V, msg);                                                                                \
        LAGraph_Free((void **)&SI, msg);                                                                               \
        LAGraph_Free((void **)&SJ, msg);                                                                               \
        LAGraph_Free((void **)&SX, msg);                                                                               \
        LAGraph_Free((void **)&parent, msg);                                                                           \
        LAGraph_Free((void **)&partner, msg);                                                                          \
        LAGraph_Free((void **)&weight, msg);                                                                           \
        GrB_free(&f);                                                                                                  \
        GrB_free(&i);                                                                                                  \
        GrB_free(&t);                                                                                                  \
        GrB_free(&edge);                                                                                               \
        GrB_free(&cedge);                                                                                              \
        GrB_free(&mask);                                                                                               \
        GrB_free(&index);                                                                                              \
        GrB_free(&comb);                                                                                               \
        GrB_free(&combMin);                                                                                            \
        GrB_free(&fst);                                                                                                \
        GrB_free(&snd);                                                                                                \
        GrB_free(&s1);                                                                                                 \
        GrB_free(&s2);                                                                                                 \
    }

//****************************************************************************
int LAGraph_msf(GrB_Matrix *result, // output: an unsymmetrical matrix, the spanning forest
                GrB_Matrix A,       // input matrix
                bool sanitize,      // if true, ensure A is symmetric
                char *msg) {

    LG_CLEAR_MSG;

#if !LAGRAPH_SUITESPARSE
    return (GrB_NOT_IMPLEMENTED);
#else

    GrB_Info info;
    GrB_Index n;
    GrB_Matrix S = NULL, T = NULL;
    GrB_Vector f = NULL, i = NULL, t = NULL, edge = NULL, cedge = NULL, mask = NULL, index = NULL;
    GrB_Index *I = NULL, *V = NULL, *SI = NULL, *SJ = NULL, *SX = NULL;

    GrB_BinaryOp comb = NULL;
    GrB_Semiring combMin = NULL;
    GrB_UnaryOp fst = NULL, snd = NULL;

    GrB_IndexUnaryOp s1 = NULL, s2 = NULL;
    if (result == NULL || A == NULL)
        return (GrB_NULL_POINTER);

    GrB_Index ncols;
    GRB_TRY(GrB_Matrix_nrows(&n, A));
    GRB_TRY(GrB_Matrix_ncols(&ncols, A));
    if (n != ncols)
        return (GrB_DIMENSION_MISMATCH);

    if (sanitize) {
        // S = A+A'
        GRB_TRY(GrB_Matrix_new(&S, GrB_UINT64, n, n));
        GRB_TRY(GrB_eWiseAdd(S, 0, 0, GrB_PLUS_UINT64, A, A, GrB_DESC_T1));
    } else {
        // Use the input as-is, and assume it is GrB_UINT64 and symmetric
        GRB_TRY(GrB_Matrix_dup(&S, A));
    }

    GRB_TRY(GrB_Matrix_new(&T, GrB_UINT64, n, n));

    GRB_TRY(GrB_Vector_new(&t, GrB_UINT64, n));
    GRB_TRY(GrB_Vector_new(&f, GrB_UINT64, n));
    GRB_TRY(GrB_Vector_new(&i, GrB_UINT64, n));
    GRB_TRY(GrB_Vector_new(&edge, GrB_UINT64, n));
    GRB_TRY(GrB_Vector_new(&cedge, GrB_UINT64, n));
    GRB_TRY(GrB_Vector_new(&mask, GrB_BOOL, n));
    GRB_TRY(GrB_Vector_new(&index, GrB_UINT64, n));

    // temporary arrays
    LG_TRY(LAGraph_Malloc((void **)&I, n, sizeof(GrB_Index), msg));
    LG_TRY(LAGraph_Malloc((void **)&V, n, sizeof(GrB_Index), msg));
    LG_TRY(LAGraph_Malloc((void **)&SI, 2 * n, sizeof(GrB_Index), msg));
    LG_TRY(LAGraph_Malloc((void **)&SJ, 2 * n, sizeof(GrB_Index), msg));
    LG_TRY(LAGraph_Malloc((void **)&SX, 2 * n, sizeof(GrB_Index), msg));

    // global arrays
    LG_TRY(LAGraph_Malloc((void **)&parent, n, sizeof(GrB_Index), msg));
    LG_TRY(LAGraph_Malloc((void **)&weight, n, sizeof(GrB_Index), msg));
    LG_TRY(LAGraph_Malloc((void **)&partner, n, sizeof(GrB_Index), msg));

    // prepare vectors
    for (GrB_Index i = 0; i < n; i++)
        I[i] = parent[i] = i;
    GRB_TRY(GrB_Vector_build(f, I, parent, n, GrB_PLUS_UINT64));
    GRB_TRY(GrB_assign(i, 0, 0, f, GrB_ALL, 0, 0));

    // semiring & monoid
    GrB_Index inf = ((uint64_t)INT_MAX << 32) ^ INT_MAX;
    GRB_TRY(GrB_BinaryOp_new(&comb, combine, GrB_UINT64, GrB_UINT64, GrB_UINT64));
    GRB_TRY(GrB_Semiring_new(&combMin, GrB_MIN_MONOID_UINT64, comb));
    GRB_TRY(GrB_UnaryOp_new(&fst, get_fst, GrB_UINT64, GrB_UINT64));
    GRB_TRY(GrB_UnaryOp_new(&snd, get_snd, GrB_UINT64, GrB_UINT64));

    // ops for GrB_select
    GrB_IndexUnaryOp_new(&s1, (void *)f1, GrB_BOOL, GrB_UINT64, GrB_UINT64);
    GrB_IndexUnaryOp_new(&s2, (void *)f2, GrB_BOOL, GrB_UINT64, GrB_UINT64);

    // the main computation
    GrB_Index nvals, diff, ntuples = 0, num;
    GRB_TRY(GrB_Matrix_nvals(&nvals, S));
    for (int iters = 1; nvals > 0; iters++) {
        // every vertex points to a root vertex at the beginning
        // edge[u] = u's minimum edge (weight and index are encoded together)
        GRB_TRY(GrB_assign(edge, 0, 0, inf, GrB_ALL, 0, 0));
        GRB_TRY(GrB_mxv(edge, 0, GrB_MIN_UINT64, combMin, S, f, 0));
        // cedge[u] = children's minimum edge  | if u is a root
        //          = (INT_MAX, u)             | otherwise
        GRB_TRY(GrB_assign(t, 0, 0, (uint64_t)INT_MAX, GrB_ALL, 0, 0));
        GRB_TRY(GrB_eWiseMult(cedge, 0, 0, comb, t, i, 0));
        LG_TRY(Reduce_assign(cedge, edge, parent, n, msg));
        // if (f[u] == u) f[u] := snd(cedge[u])  -- the index part of the edge
        GRB_TRY(GrB_eWiseMult(mask, 0, 0, GrB_EQ_UINT64, f, i, 0));
        GRB_TRY(GrB_apply(f, mask, GrB_SECOND_UINT64, snd, cedge, 0));
        // identify all the vertex pairs (u, v) where f[u] == v and f[v] == u
        // and then select the minimum of u, v as the new root;
        // if (f[f[i]] == i) f[i] = min(f[i], i)
        GRB_TRY(GrB_Vector_extractTuples(I, V, &n, f));
        GRB_TRY(GrB_extract(t, 0, 0, f, V, n, 0));
        GRB_TRY(GrB_eWiseMult(mask, 0, 0, GrB_EQ_UINT64, i, t, 0));
        GRB_TRY(GrB_assign(f, mask, GrB_MIN_UINT64, i, GrB_ALL, 0, 0));

        // five steps to generate the solution
        // 1. new roots (f[i] == i) revise their entries in cedge
        GRB_TRY(GrB_eWiseMult(mask, 0, 0, GrB_EQ_UINT64, i, f, 0));
        GRB_TRY(GrB_assign(cedge, mask, 0, inf, GrB_ALL, 0, 0));

        // 2. every vertex tries to know whether one of its edges is selected
        GRB_TRY(GrB_extract(t, 0, 0, cedge, parent, n, 0));
        GRB_TRY(GrB_eWiseMult(mask, 0, 0, GrB_EQ_UINT64, edge, t, 0));

        // 3. each root picks a vertex from its children to generate the solution
        GRB_TRY(GrB_assign(index, 0, 0, n, GrB_ALL, 0, 0));
        GRB_TRY(GrB_assign(index, mask, 0, i, GrB_ALL, 0, 0));
        GRB_TRY(GrB_assign(t, 0, 0, n, GrB_ALL, 0, 0));
        LG_TRY(Reduce_assign(t, index, parent, n, msg));
        GRB_TRY(GrB_extract(index, 0, 0, t, parent, n, 0));
        GRB_TRY(GrB_eWiseMult(mask, 0, 0, GrB_EQ_UINT64, i, index, 0));

        // 4. generate the select function (set the global pointers)
        GRB_TRY(GrB_assign(t, 0, 0, inf, GrB_ALL, 0, 0));
        GRB_TRY(GrB_apply(t, mask, 0, fst, edge, 0));
        GRB_TRY(GrB_Vector_extractTuples(I, weight, &n, t));
        GRB_TRY(GrB_assign(t, 0, 0, inf, GrB_ALL, 0, 0));
        GRB_TRY(GrB_apply(t, mask, 0, snd, edge, 0));
        GRB_TRY(GrB_Vector_extractTuples(I, partner, &n, t));
        GRB_TRY(GrB_select(T, 0, 0, s1, S, 0, 0));
        GRB_TRY(GrB_Vector_clear(t));

        // 5. the generated matrix may still have redundant edges
        //    remove the duplicates by GrB_mxv() and store them as tuples
        GRB_TRY(GrB_Vector_clear(edge));
        GRB_TRY(GrB_mxv(edge, mask, GrB_MIN_UINT64, combMin, T, i, 0));
        GRB_TRY(GrB_Vector_nvals(&num, edge));
        GRB_TRY(GrB_apply(t, 0, 0, snd, edge, 0));
        GRB_TRY(GrB_Vector_extractTuples(SI + ntuples, SJ + ntuples, &num, t));
        GRB_TRY(GrB_apply(t, 0, 0, fst, edge, 0));
        GRB_TRY(GrB_Vector_extractTuples(SI + ntuples, SX + ntuples, &num, t));
        GRB_TRY(GrB_Vector_clear(t));
        ntuples += num;

        // path halving until every vertex points on a root
        do {
            GRB_TRY(GrB_Vector_extractTuples(I, V, &n, f));
            GRB_TRY(GrB_extract(t, 0, 0, f, V, n, 0));
            GRB_TRY(GrB_eWiseMult(mask, 0, 0, GrB_NE_UINT64, f, t, 0));
            GRB_TRY(GrB_assign(f, 0, 0, t, GrB_ALL, 0, 0));
            GRB_TRY(GrB_reduce(&diff, 0, GrB_PLUS_MONOID_UINT64, mask, 0));
        } while (diff != 0);

        // remove the edges in the same connected component
        GRB_TRY(GrB_Vector_extractTuples(I, parent, &n, f));
        GRB_TRY(GrB_select(S, 0, 0, s2, S, 0, 0));
        GrB_Matrix_nvals(&nvals, S);
        if (nvals == 0)
            break;
    }
    GRB_TRY(GrB_Matrix_clear(T));
    GRB_TRY(GrB_Matrix_build(T, SI, SJ, SX, ntuples, GrB_SECOND_UINT64));
    *result = T;
    T = NULL;

    LG_FREE_ALL;
    return GrB_SUCCESS;
#endif
}