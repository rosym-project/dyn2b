#include <dyn2b/functions/linear_algebra.h>

#include <assert.h>


void la_dscal_o(
        int n,
        double alpha,
        const double *x, int incx,
        double *y, int incy)
{
    assert(x);

    for (int i = 0; i < n; i++) {
        y[i * incy] = alpha * x[i * incx];
    }
}


void la_dscal_i(
        int n,
        double alpha,
        double *x, int incx)
{
    assert(x);

    for (int i = 0; i < n; i++) {
        x[i * incx] = alpha * x[i * incx];
    }
}


void la_dgeadd_os(
        int m, int n,
        const double *a, int lda,
        const double *b, int ldb,
        double *c, int ldc)
{
    assert(a);
    assert(b);
    assert(c);
    assert(lda >= 1 && lda >= n);
    assert(ldb >= 1 && ldb >= n);
    assert(ldc >= 1 && ldc >= n);

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            c[i * ldc + j] = a[i * lda + j] + b[i * ldb + j];
        }
    }
}


void la_dgeadd_is(
        int m, int n,
        const double *a, int lda,
        double *b, int ldb)
{
    assert(a);
    assert(b);
    assert(lda >= 1 && lda >= n);
    assert(ldb >= 1 && ldb >= n);

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            b[i * ldb + j] += a[i * lda + j];
        }
    }
}


void la_daxpy_oe(
        int n,
        double alpha,
        const double *x, int incx,
        const double *y, int incy,
        double *z, int incz)
{
    assert(x);
    assert(y);
    assert(z);

    for (int i = 0; i < n; i++) {
        z[i * incz] = alpha * x[i * incx] + y[i * incy];
    }
}


void la_daxpy_ie(
        int n,
        double alpha,
        const double *x, int incx,
        double *y, int incy)
{
    assert(x);
    assert(y);
    assert(x != y);

    for (int i = 0; i < n; i++) {
        y[i * incy] = alpha * x[i * incx] + y[i * incy];
    }
}


void la_ddot(
        int n,
        double *x, int incx,
        double *y, int incy,
        double *alpha)
{
    assert(x);
    assert(y);

    *alpha = 0.0;
    for (int i = 0; i < n; i++) {
        *alpha += x[i * incx] * y[i * incy];
    }
}


void la_dcross_o(
        const double *x, int incx,
        const double *y, int incy,
        double *z, int incz)
{
    assert(x);
    assert(y);
    assert(z);
    assert(x != z);
    assert(y != z);

    // z1 = a2 b3 - a3 b2
    z[0 * incz] = x[1 * incx] * y[2 * incy] - x[2 * incx] * y[1 * incy];

    // z2 = a3 b1 - a1 b3
    z[1 * incz] = x[2 * incx] * y[0 * incy] - x[0 * incx] * y[2 * incy];

    // z3 = a1 b2 - a2 b1
    z[2 * incz] = x[0 * incx] * y[1 * incy] - x[1 * incx] * y[0 * incy];
}


void la_dcrossop(
        const double *x, int incx,
        double *a, int lda)
{
    assert(x);
    assert(a);

    a[0 * lda + 0] =          0.0; a[0 * lda + 1] = -x[2 * incx]; a[0 * lda + 2] =  x[1 * incx];
    a[1 * lda + 0] =  x[2 * incx]; a[1 * lda + 1] =          0.0; a[1 * lda + 2] = -x[0 * incx];
    a[2 * lda + 0] = -x[1 * incx]; a[2 * lda + 1] =  x[0 * incx]; a[2 * lda + 2] =          0.0;
}


void la_dgemv_nos(
        int m, int n,
        const double *a, int lda,
        const double *x, int incx,
        double *y, int incy)
{
    assert(a);
    assert(x);
    assert(y);
    assert(x != y);
    assert(lda >= 1 && lda >= n);
    assert(incx > 0);
    assert(incy > 0);

    for (int i = 0; i < n; i++) {
        double yi = 0.0;
        for (int j = 0; j < m; j++) {
            yi += a[i * lda + j] * x[j * incx];
        }
        y[i * incy] = yi;
    }
}


void la_dgemv_tos(
        int m, int n,
        const double *a, int lda,
        const double *x, int incx,
        double *y, int incy)
{
    assert(a);
    assert(x);
    assert(y);
    assert(x != y);
    assert(lda >= 1 && lda >= n);
    assert(incx > 0);
    assert(incy > 0);

    for (int i = 0; i < n; i++) {
        double yi = 0.0;
        for (int j = 0; j < m; j++) {
            yi += a[j * lda + i] * x[j * incx];
        }
        y[i * incy] = yi;
    }
}


void la_dgemv_noe(
        int m, int n,
        double alpha,
        const double *a, int lda,
        const double *x, int incx,
        double beta,
        const double *y, int incy,
        double *z, int incz)
{
    assert(a);
    assert(x);
    assert(y);
    assert(z);
    assert(x != z);
    assert(x != y);
    assert(lda >= 1 && lda >= n);
    assert(incx > 0);
    assert(incy > 0);
    assert(incz > 0);

    for (int i = 0; i < n; i++) {
        double zi = beta * y[i * incy];
        for (int j = 0; j < m; j++) {
            zi += alpha * a[i * lda + j] * x[j * incx];
        }
        z[i * incz] = zi;
    }
}


void la_dgemv_toe(
        int m, int n,
        double alpha,
        const double *a, int lda,
        const double *x, int incx,
        double beta,
        const double *y, int incy,
        double *z, int incz)
{
    assert(a);
    assert(x);
    assert(y);
    assert(z);
    assert(x != z);
    assert(x != y);
    assert(lda >= 1 && lda >= n);
    assert(incx > 0);
    assert(incy > 0);
    assert(incz > 0);

    for (int i = 0; i < n; i++) {
        double zi = beta * y[i * incy];
        for (int j = 0; j < m; j++) {
            zi += alpha * a[j * lda + i] * x[j * incx];
        }
        z[i * incz] = zi;
    }
}


void la_dgemm_nnos(
        int m, int n, int k,
        const double *a, int lda,
        const double *b, int ldb,
        double *c, int ldc)
{
    assert(a);
    assert(b);
    assert(c);
    assert(a != c);
    assert(b != c);
    assert(m >= 0);
    assert(n >= 0);
    assert(k >= 0);
    assert(lda >= 1 && lda >= k);
    assert(ldb >= 1 && ldb >= n);
    assert(ldc >= 1 && ldc >= n);

    for (int i_ = 0; i_ < m; i_++) {
        for (int j_ = 0; j_ < n; j_++) {
            double cij = 0.0;
            for (int k_ = 0; k_ < k; k_++) {
                cij += a[i_ * lda + k_] * b[k_ * ldb + j_];
            }
            c[i_ * ldc + j_] = cij;
        }
    }
}


void la_dgemm_ntos(
        int m, int n, int k,
        const double *a, int lda,
        const double *b, int ldb,
        double *c, int ldc)
{
    assert(a);
    assert(b);
    assert(c);
    assert(a != c);
    assert(b != c);
    assert(m >= 0);
    assert(n >= 0);
    assert(k >= 0);
    assert(lda >= 1 && lda >= k);
    assert(ldb >= 1 && ldb >= n);
    assert(ldc >= 1 && ldc >= n);

    for (int i_ = 0; i_ < m; i_++) {
        for (int j_ = 0; j_ < n; j_++) {
            double cij = 0.0;
            for (int k_ = 0; k_ < k; k_++) {
                cij += a[i_ * lda + k_] * b[j_ * ldb + k_];
            }
            c[i_ * ldc + j_] = cij;
        }
    }
}


void la_dgemm_tnos(
        int m, int n, int k,
        const double *a, int lda,
        const double *b, int ldb,
        double *c, int ldc)
{
    assert(a);
    assert(b);
    assert(c);
    assert(a != c);
    assert(b != c);
    assert(m >= 0);
    assert(n >= 0);
    assert(k >= 0);
    assert(lda >= 1 && lda >= m);
    assert(ldb >= 1 && ldb >= n);
    assert(ldc >= 1 && ldc >= n);

    for (int i_ = 0; i_ < m; i_++) {
        for (int j_ = 0; j_ < n; j_++) {
            double cij = 0.0;
            for (int k_ = 0; k_ < k; k_++) {
                cij += a[k_ * lda + i_] * b[k_ * ldb + j_];
            }
            c[i_ * ldc + j_] = cij;
        }
    }
}


void la_dgemm_ttos(
        int m, int n, int k,
        const double *a, int lda,
        const double *b, int ldb,
        double *c, int ldc)
{
    assert(a);
    assert(b);
    assert(c);
    assert(a != c);
    assert(b != c);
    assert(m >= 0);
    assert(n >= 0);
    assert(k >= 0);
    assert(lda >= 1 && lda >= m);
    assert(ldb >= 1 && ldb >= n);
    assert(ldc >= 1 && ldc >= n);

    for (int i_ = 0; i_ < m; i_++) {
        for (int j_ = 0; j_ < n; j_++) {
            double cij = 0.0;
            for (int k_ = 0; k_ < k; k_++) {
                cij += a[k_ * lda + i_] * b[j_ * ldb + k_];
            }
            c[i_ * ldc + j_] = cij;
        }
    }
}


void la_dgemm_nnoe(
        int m, int n, int k,
        double alpha,
        const double *a, int lda,
        const double *b, int ldb,
        double beta,
        const double *c, int ldc,
        double *d, int ldd)
{
    assert(a);
    assert(b);
    assert(c);
    assert(a != d);
    assert(b != d);
    assert(c != d);
    assert(m >= 0);
    assert(n >= 0);
    assert(k >= 0);
    assert(lda >= 1 && lda >= m);
    assert(ldb >= 1 && ldb >= n);
    assert(ldc >= 1 && ldc >= n);
    assert(ldd >= 1 && ldd >= n);

    for (int i_ = 0; i_ < m; i_++) {
        for (int j_ = 0; j_ < n; j_++) {
            double dij = beta * c[i_ * ldc + j_];
            for (int k_ = 0; k_ < k; k_++) {
                dij += alpha * a[i_ * lda + k_] * b[k_ * ldb + j_];
            }
            d[i_ * ldc + j_] = dij;
        }
    }
}


void la_dgemm_tnoe(
        int m, int n, int k,
        double alpha,
        const double *a, int lda,
        const double *b, int ldb,
        double beta,
        const double *c, int ldc,
        double *d, int ldd)
{
    assert(a);
    assert(b);
    assert(c);
    assert(a != d);
    assert(b != d);
    assert(c != d);
    assert(m >= 0);
    assert(n >= 0);
    assert(k >= 0);
    assert(lda >= 1 && lda >= m);
    assert(ldb >= 1 && ldb >= n);
    assert(ldc >= 1 && ldc >= n);
    assert(ldd >= 1 && ldd >= n);

    for (int i_ = 0; i_ < m; i_++) {
        for (int j_ = 0; j_ < n; j_++) {
            double dij = beta * c[i_ * ldc + j_];
            for (int k_ = 0; k_ < k; k_++) {
                dij += alpha * a[k_ * lda + i_] * b[k_ * ldb + j_];
            }
            d[i_ * ldc + j_] = dij;
        }
    }
}


void la_dgemm_ntoe(
        int m, int n, int k,
        double alpha,
        const double *a, int lda,
        const double *b, int ldb,
        double beta,
        const double *c, int ldc,
        double *d, int ldd)
{
    assert(a);
    assert(b);
    assert(c);
    assert(a != d);
    assert(b != d);
    assert(c != d);
    assert(m >= 0);
    assert(n >= 0);
    assert(k >= 0);
    assert(lda >= 1 && lda >= m);
    assert(ldb >= 1 && ldb >= n);
    assert(ldc >= 1 && ldc >= n);
    assert(ldd >= 1 && ldd >= n);

    for (int i_ = 0; i_ < m; i_++) {
        for (int j_ = 0; j_ < n; j_++) {
            double dij = beta * c[i_ * ldc + j_];
            for (int k_ = 0; k_ < k; k_++) {
                dij += alpha * a[i_ * lda + k_] * b[j_ * ldb + k_];
            }
            d[i_ * ldc + j_] = dij;
        }
    }
}