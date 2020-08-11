#include <dyn2b/functions/linear_algebra.h>
#include <dyn2b/types/linear_algebra.h>
#include <check.h>
#include <math.h>


#define ck_assert_flt_eq(X, Y) ck_assert_double_eq_tol(X, Y, 0.0001)


static struct matrix3x3 mat_a = {
    .row_x = { 1.0, 2.0, 3.0 },
    .row_y = { 4.0, 5.0, 6.0 },
    .row_z = { 7.0, 8.0, 9.0 } };
static struct matrix3x3 mat_b = {
    .row_x = { 2.0, 3.0,  4.0 },
    .row_y = { 5.0, 6.0,  7.0 },
    .row_z = { 8.0, 9.0, 10.0 } };


START_TEST(test_la_dscal_o)
{
    struct vector3 a = { 1.0, 2.0, 3.0 };
    struct vector3 r;

    la_dscal_o(3, -2.0, (double *)&a, 1, (double *)&r, 1);
    ck_assert_flt_eq(r.x, -2.0);
    ck_assert_flt_eq(r.y, -4.0);
    ck_assert_flt_eq(r.z, -6.0);
}
END_TEST


START_TEST(test_la_dscal_i)
{
    struct vector3 a = { 1.0, 2.0, 3.0 };

    la_dscal_i(3, -2.0, (double *)&a, 1);
    ck_assert_flt_eq(a.x, -2.0);
    ck_assert_flt_eq(a.y, -4.0);
    ck_assert_flt_eq(a.z, -6.0);
}
END_TEST


START_TEST(test_la_dgeadd_os)
{
    struct matrix3x3 c = {
        .row_x = { 1.0, 2.0, 3.0 },
        .row_y = { 4.0, 5.0, 6.0 },
        .row_z = { 7.0, 8.0, 9.0 } };

    struct matrix3x3 r1 = {
        .row_x = {  3.0,  5.0,  7.0 },
        .row_y = {  9.0, 11.0, 13.0 },
        .row_z = { 15.0, 17.0, 19.0 } };


    // Operate in in-place mode
    la_dgeadd_os(3, 3,
            (double *)&c, 3,
            (double *)&mat_b, 3,
            (double *)&c, 3);
    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 1; j++) {
            ck_assert_flt_eq(c.row[i].data[j], r1.row[i].data[j]);
        }
    }


    // Operate in out-of-place mode
    la_dgeadd_os(3, 3,
            (double *)&mat_a, 3,
            (double *)&mat_b, 3,
            (double *)&c, 3);
    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 1; j++) {
            ck_assert_flt_eq(c.row[i].data[j], r1.row[i].data[j]);
        }
    }

    struct matrix3x3 r2 = {
        .row_x = {  2.0,  4.0,  6.0 },
        .row_y = {  8.0, 10.0, 12.0 },
        .row_z = { 14.0, 16.0, 18.0 } };

    la_dgeadd_os(3, 3,
            (double *)&mat_a, 3,
            (double *)&mat_a, 3,
            (double *)&c, 3);
    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 1; j++) {
            ck_assert_flt_eq(c.row[i].data[j], r2.row[i].data[j]);
        }
    }
}
END_TEST


START_TEST(test_la_dgeadd_is)
{
    struct matrix3x3 c = {
        .row_x = { 1.0, 2.0, 3.0 },
        .row_y = { 4.0, 5.0, 6.0 },
        .row_z = { 7.0, 8.0, 9.0 } };

    struct matrix3x3 r1 = {
        .row_x = {  3.0,  5.0,  7.0 },
        .row_y = {  9.0, 11.0, 13.0 },
        .row_z = { 15.0, 17.0, 19.0 } };


    la_dgeadd_is(3, 3,
            (double *)&mat_b, 3,
            (double *)&c, 3);
    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 1; j++) {
            ck_assert_flt_eq(c.row[i].data[j], r1.row[i].data[j]);
        }
    }
}
END_TEST


START_TEST(test_la_daxpy_oe)
{
    struct vector3 a = { 1.0, 2.0, 3.0 };
    struct vector3 b = { 2.0, 3.0, 4.0 };
    struct vector3 r;

    la_daxpy_oe(3,
            -2.0, (double *)&a, 1,
            (double *)&b, 1,
            (double *)&r, 1);
    ck_assert_flt_eq(r.x,  0.0);
    ck_assert_flt_eq(r.y, -1.0);
    ck_assert_flt_eq(r.z, -2.0);
}
END_TEST


START_TEST(test_la_daxpy_ie)
{
    struct vector3 a = { 1.0, 2.0, 3.0 };
    struct vector3 b = { 2.0, 3.0, 4.0 };

    la_daxpy_ie(3,
            -2.0, (double *)&a, 1,
            (double *)&b, 1);
    ck_assert_flt_eq(b.x,  0.0);
    ck_assert_flt_eq(b.y, -1.0);
    ck_assert_flt_eq(b.z, -2.0);
}
END_TEST


START_TEST(test_la_ddot)
{
    struct vector3 a = { 1.0, 2.0, 3.0 };
    struct vector3 b = { 2.0, 3.0, 4.0 };
    double r;

    la_ddot(3,
            (double *)&a, 1,
            (double *)&b, 1,
            &r);
    ck_assert_flt_eq(r, 20.0);
}
END_TEST


START_TEST(test_la_dcross_o)
{
    struct vector3 a = { 1.0, 2.0, 3.0 };
    struct vector3 b = { 2.0, 3.0, 4.0 };
    struct vector3 r;

    la_dcross_o((double *)&a, 1, (double *)&a, 1, (double *)&r, 1);
    ck_assert_flt_eq(r.x, 0.0);
    ck_assert_flt_eq(r.y, 0.0);
    ck_assert_flt_eq(r.z, 0.0);

    la_dcross_o((double *)&a, 1, (double *)&b, 1, (double *)&r, 1);
    ck_assert_flt_eq(r.x, -1.0);
    ck_assert_flt_eq(r.y,  2.0);
    ck_assert_flt_eq(r.z, -1.0);
}
END_TEST


START_TEST(test_la_dcrossop)
{
    struct vector3 a = { 1.0, 2.0, 3.0 };
    struct matrix3x3 b;

    struct matrix3x3 r = {
        .row_x = {  0.0, -3.0,  2.0 },
        .row_y = {  3.0,  0.0, -1.0 },
        .row_z = { -2.0,  1.0,  0.0 } };

    la_dcrossop((double *)&a, 1, (double *)&b, 3);
    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 1; j++) {
            ck_assert_flt_eq(b.row[i].data[j], r.row[i].data[j]);
        }
    }
}
END_TEST


START_TEST(test_la_dgemv_nos)
{
    struct vector3 b = { 2.0, 3.0, 4.0 };
    struct vector3 r;


    la_dgemv_nos(3, 3, (double *)&mat_a, 3, (double *)&b, 1, (double *)&r, 1);
    ck_assert_flt_eq(r.x, 20.0);
    ck_assert_flt_eq(r.y, 47.0);
    ck_assert_flt_eq(r.z, 74.0);


    double arr_a[2][3] = { { 1.0, 2.0, 3.0 },
                           { 4.0, 5.0, 6.0 } };
    double arr_b[3] = { 2.0, 3.0, 4.0 };
    double arr_r[2];

    la_dgemv_nos(3, 2, &arr_a[0][0], 3, &arr_b[0], 1, &arr_r[0], 1);
    ck_assert_flt_eq(arr_r[0], 20.0);
    ck_assert_flt_eq(arr_r[1], 47.0);
}
END_TEST


START_TEST(test_la_dgemv_tos)
{
    struct vector3 b = { 2.0, 3.0, 4.0 };
    struct vector3 r;


    la_dgemv_tos(3, 3, (double *)&mat_a, 3, (double *)&b, 1, (double *)&r, 1);
    ck_assert_flt_eq(r.x, 42.0);
    ck_assert_flt_eq(r.y, 51.0);
    ck_assert_flt_eq(r.z, 60.0);


    double arr_a[3][2] = { { 1.0, 4.0 },
                           { 2.0, 5.0 },
                           { 3.0, 6.0 } };
    double arr_b[3] = { 2.0, 3.0, 4.0 };
    double arr_r[2];

    la_dgemv_tos(3, 2, &arr_a[0][0], 2, &arr_b[0], 1, &arr_r[0], 1);
    ck_assert_flt_eq(arr_r[0], 20.0);
    ck_assert_flt_eq(arr_r[1], 47.0);
}
END_TEST


START_TEST(test_la_dgemv_noe)
{
    struct vector3 b = { 2.0, 3.0, 4.0 };
    struct vector3 c = { 3.0, 4.0, 5.0 };
    struct vector3 r;


    la_dgemv_noe(3, 3,
            2.0, (double *)&mat_a, 3, (double *)&b, 1,
            3.0, (double *)&c, 1,
            (double *)&r, 1);
    ck_assert_flt_eq(r.x,  49.0);
    ck_assert_flt_eq(r.y, 106.0);
    ck_assert_flt_eq(r.z, 163.0);


    double arr_a[2][3] = { { 1.0, 2.0, 3.0 },
                           { 4.0, 5.0, 6.0 } };
    double arr_b[3] = { 2.0, 3.0, 4.0 };
    double arr_c[2] = { 3.0, 4.0 };
    double arr_r[2];

    la_dgemv_noe(3, 2,
            2.0, &arr_a[0][0], 3, &arr_b[0], 1,
            3.0, &arr_c[0], 1,
            &arr_r[0], 1);
    ck_assert_flt_eq(arr_r[0],  49.0);
    ck_assert_flt_eq(arr_r[1], 106.0);
}
END_TEST


START_TEST(test_la_dgemv_toe)
{
    struct vector3 b = { 2.0, 3.0, 4.0 };
    struct vector3 c = { 3.0, 4.0, 5.0 };
    struct vector3 r;


    la_dgemv_toe(3, 3,
            2.0, (double *)&mat_a, 3, (double *)&b, 1,
            3.0, (double *)&c, 1,
            (double *)&r, 1);
    ck_assert_flt_eq(r.x,  93.0);
    ck_assert_flt_eq(r.y, 114.0);
    ck_assert_flt_eq(r.z, 135.0);


    double arr_a[3][2] = { { 1.0, 4.0 },
                           { 2.0, 5.0 },
                           { 3.0, 6.0 } };
    double arr_b[3] = { 2.0, 3.0, 4.0 };
    double arr_c[2] = { 3.0, 4.0 };
    double arr_r[2];

    la_dgemv_toe(3, 2,
            2.0, &arr_a[0][0], 2, &arr_b[0], 1,
            3.0, &arr_c[0], 1,
            &arr_r[0], 1);
    ck_assert_flt_eq(arr_r[0],  49.0);
    ck_assert_flt_eq(arr_r[1], 106.0);
}
END_TEST


START_TEST(test_la_dgemm_nnos)
{
    struct matrix3x3 c;


    struct matrix3x3 r1 = {
        .row_x = {  30.0,  36.0,  42.0 },
        .row_y = {  66.0,  81.0,  96.0 },
        .row_z = { 102.0, 126.0, 150.0 } };

    la_dgemm_nnos(3, 3, 3, (double *)&mat_a, 3, (double *)&mat_a, 3, (double *)&c, 3);
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            ck_assert_flt_eq(c.row[i].data[j], r1.row[i].data[j]);
        }
    }


    struct matrix3x3 r2 = {
        .row_x = {  36.0,  42.0,  48.0 },
        .row_y = {  81.0,  96.0, 111.0 },
        .row_z = { 126.0, 150.0, 174.0 } };

    la_dgemm_nnos(3, 3, 3, (double *)&mat_a, 3, (double *)&mat_b, 3, (double *)&c, 3);
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            ck_assert_flt_eq(c.row[i].data[j], r2.row[i].data[j]);
        }
    }


    double arr_a[2][3] = { { 1.0, 2.0, 3.0 },
                           { 4.0, 5.0, 6.0 } };
    double arr_b[3][1] = { { 7.0 }, { 8.0 }, { 9.0 } };
    double arr_c[2][1];

    double arr_r[2][1] = { { 50.0 },
                             { 122.0 } };

    la_dgemm_nnos(2, 1, 3, &arr_a[0][0], 3, &arr_b[0][0], 1, &arr_c[0][0], 1);
    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 1; j++) {
            ck_assert_flt_eq(arr_c[i][j], arr_r[i][j]);
        }
    }
}
END_TEST


START_TEST(test_la_dgemm_ntos)
{
    struct matrix3x3 c;


    struct matrix3x3 r1 = {
        .row_x = { 14.0,  32.0,  50.0 },
        .row_y = { 32.0,  77.0, 122.0 },
        .row_z = { 50.0, 122.0, 194.0 } };

    la_dgemm_ntos(3, 3, 3, (double *)&mat_a, 3, (double *)&mat_a, 3, (double *)&c, 3);
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            ck_assert_flt_eq(c.row[i].data[j], r1.row[i].data[j]);
        }
    }


    struct matrix3x3 r2 = {
        .row_x = { 20.0,  38.0,  56.0 },
        .row_y = { 47.0,  92.0, 137.0 },
        .row_z = { 74.0, 146.0, 218.0 } };

    la_dgemm_ntos(3, 3, 3, (double *)&mat_a, 3, (double *)&mat_b, 3, (double *)&c, 3);
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            ck_assert_flt_eq(c.row[i].data[j], r2.row[i].data[j]);
        }
    }


    double arr_a[2][3] = { { 1.0, 2.0, 3.0 },
                           { 4.0, 5.0, 6.0 } };
    double arr_b[1][3] = { { 7.0, 8.0, 9.0 } };
    double arr_c[2][1];

    double arr_r[2][1] = { {  50.0 },
                           { 122.0 } };

    la_dgemm_ntos(2, 1, 3, &arr_a[0][0], 3, &arr_b[0][0], 3, &arr_c[0][0], 1);
    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 1; j++) {
            ck_assert_flt_eq(arr_c[i][j], arr_r[i][j]);
        }
    }
}
END_TEST


START_TEST(test_la_dgemm_tnos)
{
    struct matrix3x3 c;


    struct matrix3x3 r1 = {
        .row_x = { 66.0,  78.0,  90.0 },
        .row_y = { 78.0,  93.0, 108.0 },
        .row_z = { 90.0, 108.0, 126.0 } };

    la_dgemm_tnos(3, 3, 3, (double *)&mat_a, 3, (double *)&mat_a, 3, (double *)&c, 3);
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            ck_assert_flt_eq(c.row[i].data[j], r1.row[i].data[j]);
        }
    }


    struct matrix3x3 r2 = {
        .row_x = {  78.0,  90.0, 102.0 },
        .row_y = {  93.0, 108.0, 123.0 },
        .row_z = { 108.0, 126.0, 144.0 } };

    la_dgemm_tnos(3, 3, 3, (double *)&mat_a, 3, (double *)&mat_b, 3, (double *)&c, 3);
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            ck_assert_flt_eq(c.row[i].data[j], r2.row[i].data[j]);
        }
    }


    double arr_a[3][2] = { { 1.0, 4.0 },
                           { 2.0, 5.0 },
                           { 3.0, 6.0 }};
    double arr_b[3][1] = { { 7.0 },
                           { 8.0 },
                           { 9.0 } };
    double arr_c[2][1];

    double arr_r[2][1] = { {  50.0 },
                           { 122.0 } };

    la_dgemm_tnos(2, 1, 3, &arr_a[0][0], 2, &arr_b[0][0], 1, &arr_c[0][0], 1);
    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 1; j++) {
            ck_assert_flt_eq(arr_c[i][j], arr_r[i][j]);
        }
    }
}
END_TEST


START_TEST(test_la_dgemm_ttos)
{
    struct matrix3x3 c;


    struct matrix3x3 r1 = {
        .row_x = { 30.0, 66.0, 102.0 },
        .row_y = { 36.0, 81.0, 126.0 },
        .row_z = { 42.0, 96.0, 150.0 } };

    la_dgemm_ttos(3, 3, 3, (double *)&mat_a, 3, (double *)&mat_a, 3, (double *)&c, 3);
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            ck_assert_flt_eq(c.row[i].data[j], r1.row[i].data[j]);
        }
    }


    struct matrix3x3 r2 = {
        .row_x = { 42.0,  78.0, 114.0 },
        .row_y = { 51.0,  96.0, 141.0 },
        .row_z = { 60.0, 114.0, 168.0 } };

    la_dgemm_ttos(3, 3, 3, (double *)&mat_a, 3, (double *)&mat_b, 3, (double *)&c, 3);
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            ck_assert_flt_eq(c.row[i].data[j], r2.row[i].data[j]);
        }
    }


    double arr_a[3][2] = { { 1.0, 4.0 },
                           { 2.0, 5.0 },
                           { 3.0, 6.0 } };
    double arr_b[1][3] = { { 7.0, 8.0, 9.0 } };
    double arr_c[2][1];

    double arr_r[2][1] = { {  50.0 },
                           { 122.0 } };

    la_dgemm_ttos(2, 1, 3, &arr_a[0][0], 2, &arr_b[0][0], 3, &arr_c[0][0], 1);
    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 1; j++) {
            ck_assert_flt_eq(arr_c[i][j], arr_r[i][j]);
        }
    }
}
END_TEST


START_TEST(test_la_dgemm_nnoe)
{
    struct matrix3x3 c = {
        .row_x = { 3.0, 4.0, 5.0 },
        .row_y = { 4.0, 5.0, 6.0 },
        .row_z = { 6.0, 7.0, 8.0 } };
    struct matrix3x3 d;


    struct matrix3x3 r = {
        .row_x = {  81.0,  96.0, 111.0 },
        .row_y = { 174.0, 207.0, 240.0 },
        .row_z = { 270.0, 321.0, 372.0 } };

    la_dgemm_nnoe(3, 3, 3,
            2.0, (double *)&mat_a, 3, (double *)&mat_b, 3,
            3.0, (double *)&c, 3,
            (double *)&d, 3);
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            ck_assert_flt_eq(d.row[i].data[j], r.row[i].data[j]);
        }
    }


    double arr_a[2][3] = { { 1.0, 2.0, 3.0 },
                           { 4.0, 5.0, 6.0 } };
    double arr_b[3][1] = { { 7.0 }, { 8.0 }, { 9.0 } };
    double arr_c[2][1] = { { 2.0 },
                           { 3.0 } };
    double arr_d[2][1];

    double arr_r[2][1] = { { 106.0 },
                           { 253.0 } };

    la_dgemm_nnoe(2, 1, 3,
            2.0, &arr_a[0][0], 3, &arr_b[0][0], 1,
            3.0, &arr_c[0][0], 1,
            &arr_d[0][0], 1);
    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 1; j++) {
            ck_assert_flt_eq(arr_d[i][j], arr_r[i][j]);
        }
    }
}
END_TEST


START_TEST(test_la_dgemm_tnoe)
{
    struct matrix3x3 c = {
        .row_x = { 3.0, 4.0, 5.0 },
        .row_y = { 4.0, 5.0, 6.0 },
        .row_z = { 6.0, 7.0, 8.0 } };
    struct matrix3x3 d;


    struct matrix3x3 r = {
        .row_x = { 165.0, 192.0, 219.0 },
        .row_y = { 198.0, 231.0, 264.0 },
        .row_z = { 234.0, 273.0, 312.0 } };

    la_dgemm_tnoe(3, 3, 3,
            2.0, (double *)&mat_a, 3, (double *)&mat_b, 3,
            3.0, (double *)&c, 3,
            (double *)&d, 3);
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            ck_assert_flt_eq(d.row[i].data[j], r.row[i].data[j]);
        }
    }


    double arr_a[3][2] = { { 1.0, 4.0 },
                           { 2.0, 5.0 },
                           { 3.0, 6.0 } };
    double arr_b[3][1] = { { 7.0 }, { 8.0 }, { 9.0 } };
    double arr_c[2][1] = { { 2.0 },
                           { 3.0 } };
    double arr_d[2][1];

    double arr_r[2][1] = { { 106.0 },
                           { 253.0 } };

    la_dgemm_tnoe(2, 1, 3,
            2.0, &arr_a[0][0], 2, &arr_b[0][0], 1,
            3.0, &arr_c[0][0], 1,
            &arr_d[0][0], 1);
    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 1; j++) {
            ck_assert_flt_eq(arr_d[i][j], arr_r[i][j]);
        }
    }
}
END_TEST


START_TEST(test_la_dgemm_ntoe)
{
    struct matrix3x3 c = {
        .row_x = { 3.0, 4.0, 5.0 },
        .row_y = { 4.0, 5.0, 6.0 },
        .row_z = { 6.0, 7.0, 8.0 } };
    struct matrix3x3 d;


    struct matrix3x3 r = {
        .row_x = {  49.0,  88.0, 127.0 },
        .row_y = { 106.0, 199.0, 292.0 },
        .row_z = { 166.0, 313.0, 460.0 } };

    la_dgemm_ntoe(3, 3, 3,
            2.0, (double *)&mat_a, 3, (double *)&mat_b, 3,
            3.0, (double *)&c, 3,
            (double *)&d, 3);
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            ck_assert_flt_eq(d.row[i].data[j], r.row[i].data[j]);
        }
    }


    double arr_a[2][3] = { { 1.0, 2.0, 3.0 },
                           { 4.0, 5.0, 6.0 } };
    double arr_b[1][3] = { { 7.0, 8.0, 9.0 } };
    double arr_c[2][1] = { { 2.0 },
                           { 3.0 } };
    double arr_d[2][1];

    double arr_r[2][1] = { { 106.0 },
                           { 253.0 } };

    la_dgemm_ntoe(2, 1, 3,
            2.0, &arr_a[0][0], 3, &arr_b[0][0], 3,
            3.0, &arr_c[0][0], 1,
            &arr_d[0][0], 1);
    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 1; j++) {
            ck_assert_flt_eq(arr_d[i][j], arr_r[i][j]);
        }
    }
}
END_TEST


TCase *linear_algebra_test()
{
    TCase *tc = tcase_create("LinearAlgebra");

    tcase_add_test(tc, test_la_dscal_o);
    tcase_add_test(tc, test_la_dscal_i);
    tcase_add_test(tc, test_la_dgeadd_os);
    tcase_add_test(tc, test_la_dgeadd_is);
    tcase_add_test(tc, test_la_daxpy_oe);
    tcase_add_test(tc, test_la_daxpy_ie);
    tcase_add_test(tc, test_la_ddot);
    tcase_add_test(tc, test_la_dcross_o);
    tcase_add_test(tc, test_la_dcrossop);
    tcase_add_test(tc, test_la_dgemv_nos);
    tcase_add_test(tc, test_la_dgemv_tos);
    tcase_add_test(tc, test_la_dgemv_noe);
    tcase_add_test(tc, test_la_dgemv_toe);
    tcase_add_test(tc, test_la_dgemm_nnos);
    tcase_add_test(tc, test_la_dgemm_ntos);
    tcase_add_test(tc, test_la_dgemm_tnos);
    tcase_add_test(tc, test_la_dgemm_ttos);
    tcase_add_test(tc, test_la_dgemm_nnoe);
    tcase_add_test(tc, test_la_dgemm_tnoe);
    tcase_add_test(tc, test_la_dgemm_ntoe);

    return tc;
}
