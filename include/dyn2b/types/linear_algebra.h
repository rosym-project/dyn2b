#ifndef DYN2B_TYPES_LINEAR_ALGEBRA_H
#define DYN2B_TYPES_LINEAR_ALGEBRA_H

#ifdef __cplusplus
extern "C" {
#endif


struct vector3
{
    union {
        struct {
            double x;
            double y;
            double z;
        };
        double data[3];
    };
};

struct matrix3x3
{
    union {
        struct {
            struct vector3 row_x;
            struct vector3 row_y;
            struct vector3 row_z;
        };
        struct vector3 row[3];
    };
};

#ifdef __cplusplus
}
#endif

#endif
