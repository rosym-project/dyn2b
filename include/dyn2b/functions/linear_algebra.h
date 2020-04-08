#ifndef DYN2B_FUNCTIONS_LINEAR_ALGEBRA_H
#define DYN2B_FUNCTIONS_LINEAR_ALGEBRA_H

#ifdef __cplusplus
extern "C" {
#endif

/**
 * BLAS-, LAPACK-, MAGMA-inspired linear algebra operations
 * 
 * The following naming conventions apply, inspired by other BLAS
 * implementations:
 * la_{precision}{matrix_type}{operation}{*flags}
 * 
 * The flags include:
 * - s/e: _s_imple = reduced operation w.r.t. BLAS vs _e_xtended = BLAS operations
 * - i/o: _i_n-place vs. _o_ut-of-place
 *   https://github.com/xianyi/OpenBLAS/wiki/OpenBLAS-Extensions
 * - t/n: transposed vs. non-transposed
 *   https://blasfeo.syscop.de/docs/naming/#additional-flags
 * - l/u: lower vs. upper
 *   https://blasfeo.syscop.de/docs/naming/#additional-flags
 */


/**
 * Vector scaling (out-of-place)
 * y = alpha * x
 * 
 * alpha: scalar
 * x: m x 1
 * y: m x 1
 * 
 * Reference:
 * - https://www.netlib.org/lapack/explore-html/d4/dd0/dscal_8f.html
 */
void la_dscal_o(
        int n,
        double alpha,
        const double *x, int incx,
        double *y, int incy);

/**
 * Vector scaling (in-place)
 * x = alpha * x
 * 
 * alpha: scalar
 * x: m x 1
 * 
 * Reference:
 * - https://www.netlib.org/lapack/explore-html/d4/dd0/dscal_8f.html
 */
void la_dscal_i(
        int n,
        double alpha,
        double *x, int incx);

/**
 * General matrix addition (out-of-place, simple)
 * C = A + B
 * 
 * A: m x n
 * B: m x n
 * C: m x n
 *
 * Reference:
 * - https://icl.cs.utk.edu/projectsfiles/magma/doxygen/group__magma__geadd.html
 * - https://github.com/xianyi/OpenBLAS/blob/develop/kernel/generic/geadd.c
 * - https://github.com/flame/blis/blob/master/docs/BLISTypedAPI.md#addm
 */
void la_dgeadd_os(
        int m, int n,
        const double *a, int lda,
        const double *b, int ldb,
        double *c, int ldc);

/**
 * General matrix addition (in-place, simple)
 * B = A + B
 * 
 * A: m x n
 * B: m x n
 *
 * Reference:
 * - https://icl.cs.utk.edu/projectsfiles/magma/doxygen/group__magma__geadd.html
 * - https://github.com/xianyi/OpenBLAS/blob/develop/kernel/generic/geadd.c
 * - https://github.com/flame/blis/blob/master/docs/BLISTypedAPI.md#addm
 */
void la_dgeadd_is(
        int m, int n,
        const double *a, int lda,
        double *b, int ldb);

/**
 * Vector addition (out-of-place, extended)
 * z = alpha * x + y
 *
 * alpha: scalar
 * x: n x 1
 * y: n x 1
 * z: n x 1
 *
 * Reference:
 * - https://www.netlib.org/lapack/explore-html/d9/dcd/daxpy_8f.html
 * - https://icl.cs.utk.edu/projectsfiles/magma/doxygen/group__magma__axpy.html
 */
void la_daxpy_oe(
        int n,
        double alpha,
        const double *x, int incx,
        const double *y, int incy,
        double *z, int incz);

/**
 * Vector addition (in-place, extended)
 * y = alpha * x + y
 *
 * alpha: scalar
 * x: n x 1
 * y: n x 1
 *
 * Reference:
 * - https://www.netlib.org/lapack/explore-html/d9/dcd/daxpy_8f.html
 * - https://icl.cs.utk.edu/projectsfiles/magma/doxygen/group__magma__axpy.html
 */
void la_daxpy_ie(
        int n,
        double alpha,
        const double *x, int incx,
        double *y, int incy);

/**
 * Vector dot product.
 * alpha = x^T * y
 * 
 * x: n x 1
 * y: n x 1
 * alpha: scalar
 *
 * Reference:
 * - https://www.netlib.org/lapack/explore-html/d5/df6/ddot_8f.html
 * - https://icl.cs.utk.edu/projectsfiles/magma/doxygen/group__magma____dot.html
 * - https://github.com/xianyi/OpenBLAS/blob/develop/kernel/generic/dot.c
 */
void la_ddot(
        int n,
        double *x, int incx,
        double *y, int incy,
        double *alpha);


/**
 * Vector cross product (out-of-place).
 * z = x x y
 * 
 * a: 3 x 1
 * b: 3 x 1
 * c: 3 x 1
 */
void la_dcross_o(
        const double *x, int incx,
        const double *y, int incy,
        double *z, int incz);

/**
 * Cross operator that maps a vector to a skew-symmetric matrix.
 * a = [x]_x
 * 
 * x: 3 x 1
 * a: 3 x 3
 */
void la_dcrossop(
        const double *x, int incx,
        double *a, int lda);

/**
 * Matrix-vector product (out-of-place, simple).
 * y = A * x
 *
 * A: m x n
 * x: n x 1
 * y: m x 1
 * 
 * Reference:
 * - https://www.netlib.org/lapack/explore-html/dc/da8/dgemv_8f.html
 */
void la_dgemv_nos(
        int m, int n,
        const double *a, int lda,
        const double *x, int incx,
        double *y, int incy);

/**
 * Matrix-vector product (out-of-place, simple)
 * y = A^T * x 
 *
 * A: m x n
 * x: m x 1
 * y: n x 1
 * 
 * Reference:
 * - https://www.netlib.org/lapack/explore-html/dc/da8/dgemv_8f.html
 */
void la_dgemv_tos(
        int m, int n,
        const double *a, int lda,
        const double *x, int incx,
        double *y, int incy);

/**
 * Matrix-vector product (out-of-place, extended).
 * z = alpha * A * x + beta * y
 *
 * alpha: scalar
 * A: m x n
 * x: n x 1
 * beta: scalar
 * y: m x 1
 * z: m x 1
 * 
 * Reference:
 * - https://www.netlib.org/lapack/explore-html/dc/da8/dgemv_8f.html
 */
void la_dgemv_noe(
        int m, int n,
        double alpha,
        const double *a, int lda,
        const double *x, int incx,
        double beta,
        const double *y, int incy,
        double *z, int incz);

/**
 * Matrix-vector product (out-of-place, extended).
 * r = alpha * A^T * x + beta * y
 *
 * alpha: scalar
 * A: m x n
 * x: m x 1
 * beta: scalar
 * y: n x 1
 * z: n x 1
 * 
 * Reference:
 * - https://www.netlib.org/lapack/explore-html/dc/da8/dgemv_8f.html
 */
void la_dgemv_toe(
        int m, int n,
        double alpha,
        const double *a, int lda,
        const double *x, int incx,
        double beta,
        const double *y, int incy,
        double *z, int incz);

/**
 * Matrix-matrix product (out-of-place, simple)
 * C = A * B
 *
 * A: m x k
 * B: k x n
 * C: m x n
 * 
 * Reference:
 * - https://www.netlib.org/lapack/explore-html/d7/d2b/dgemm_8f.html
 */
void la_dgemm_nnos(
        int m, int n, int k,
        const double *a, int lda,
        const double *b, int ldb,
        double *c, int ldc);

/**
 * Matrix-matrix product (out-of-place, simple)
 * C = A * B^T
 *
 * A: m x k
 * B: n x k
 * C: m x n
 * 
 * Reference:
 * - https://www.netlib.org/lapack/explore-html/d7/d2b/dgemm_8f.html
 */
void la_dgemm_ntos(
        int m, int n, int k,
        const double *a, int lda,
        const double *b, int ldb,
        double *c, int ldc);

/**
 * Matrix-matrix product (out-of-place, simple)
 * C = A^T * B
 *
 * A: k x m
 * B: k x n
 * C: m x n
 * 
 * Reference:
 * - https://www.netlib.org/lapack/explore-html/d7/d2b/dgemm_8f.html
 */
void la_dgemm_tnos(
        int m, int n, int k,
        const double *a, int lda,
        const double *b, int ldb,
        double *c, int ldc);

/**
 * Matrix-matrix product (out-of-place, simple)
 * C = A^T * B^T
 *
 * A: k x m
 * B: n x k
 * C: m x n
 * 
 * Reference:
 * - https://www.netlib.org/lapack/explore-html/d7/d2b/dgemm_8f.html
 */
void la_dgemm_ttos(
        int m, int n, int k,
        const double *a, int lda,
        const double *b, int ldb,
        double *c, int ldc);

/**
 * Matrix-vector product (out-of-place, extended).
 * D = alpha * A * B + beta * C
 *
 * alpha: scalar
 * A: m x k
 * B: k x n
 * beta: scalar
 * C: m x n
 * D: m x n
 * 
 * Reference:
 * - https://www.netlib.org/lapack/explore-html/d7/d2b/dgemm_8f.html
 */
void la_dgemm_nnoe(
        int m, int n, int k,
        double alpha,
        const double *a, int lda,
        const double *b, int ldb,
        double beta,
        const double *c, int ldc,
        double *d, int ldd);

/**
 * Matrix-vector product (out-of-place, extended).
 * D = alpha * A^T * B + beta * C
 *
 * alpha: scalar
 * A: m x k
 * B: k x n
 * beta: scalar
 * C: m x n
 * D: m x n
 * 
 * Reference:
 * - https://www.netlib.org/lapack/explore-html/d7/d2b/dgemm_8f.html
 */
void la_dgemm_tnoe(
        int m, int n, int k,
        double alpha,
        const double *a, int lda,
        const double *b, int ldb,
        double beta,
        const double *c, int ldc,
        double *d, int ldd);

/**
 * Matrix-vector product (out-of-place, extended).
 * D = alpha * A * B^T + beta * C
 *
 * alpha: scalar
 * A: m x k
 * B: n x k
 * beta: scalar
 * C: m x n
 * D: m x n
 * 
 * Reference:
 * - https://www.netlib.org/lapack/explore-html/d7/d2b/dgemm_8f.html
 */
void la_dgemm_ntoe(
        int m, int n, int k,
        double alpha,
        const double *a, int lda,
        const double *b, int ldb,
        double beta,
        const double *c, int ldc,
        double *d, int ldd);



/**
 * Rank-one update of a symmetric matrix (lower triangular, out-of-place)
 * B = alpha * x * x^T + A
 * 
 * A: n x n
 * alpha: scalar
 * x: n x 1
 * B: n x n
 * 
 * Reference:
 * - https://www.netlib.org/lapack/explore-html/d3/d60/dsyr_8f.html
 * - https://icl.cs.utk.edu/projectsfiles/magma/doxygen/group__magma__syr.html
 * - https://github.com/xianyi/OpenBLAS/blob/develop/interface/syr.c
 */
//void la_dsyr_lo(
//        int n,
//        double alpha,
//        const double *x, int incx,
//        const double *a, int lda,
//        double *b, int ldb);

/**
 * LDL^T factorization of a symmetric matrix (out-of-place)
 * A = L * D * L^T
 * 
 * The factors (L, D) are stored in B.
 * 
 * Reference:
 * - https://www.netlib.org/lapack/explore-html/dd/df4/dsytrf_8f.html
 * - https://www.netlib.org/lapack/explore-html/d0/d9a/dsytrs_8f.html
 * - https://blasfeo.syscop.de/docs/naming/#factorization
 */
//void la_dsytrf_los(
//    int n,
//    const double *a, int lda,
//    double *b, int ldb);

/**
 * Rank-one update of an LDL^T factorization of a symmetric matrix (out-of-place)
 * (LDL^T)_B = alpha * x * x^T + (LDL^T)_A
 *
 * A: n x n
 * alpha: scalar
 * x: n x 1
 * B: n x n
 *
 * Reference:
 * - https://www.netlib.org/lapack/explore-html/dd/df4/dsytrf_8f.html
 * - https://www.netlib.org/lapack/explore-html/d3/d60/dsyr_8f.html
 */
void la_dsytrfr_lo(
        int n,
        double alpha,
        const double *x, int incx,
        const double *a, int lda,
        double *b, int ldb);

/**
 * Solve system of equations (lower triangular, no transpose, diagonal unit, out-of-place)
 * A * x = b <=> x = A^{-1} b
 * 
 * "Forward substitution"
 *
 * A: lda x n
 * b: n x 1
 * x: n x 1
 *
 * Reference:
 * - https://www.netlib.org/lapack/explore-html/d6/d96/dtrsv_8f.html
 */
void la_trsv_lnd(
        int n,
        const double *a, int lda,
        const double *b, int incb,
        double *x, int incx);

/**
 * Solve system of equations (lower triangular, transpose, diagonal unit, out-of-place)
 * A^T * x = b <=> x = (A^T)^{-1} b
 * 
 * "Back substitution"
 *
 * A: n x lda
 * b: n x 1
 * x: n x 1
 *
 * Reference:
 * - https://www.netlib.org/lapack/explore-html/d6/d96/dtrsv_8f.html
 */
void la_trsv_ltd(
        int n,
        const double *a, int lda,
        const double *b, int incb,
        double *x, int incx);

#ifdef __cplusplus
}
#endif

#endif
