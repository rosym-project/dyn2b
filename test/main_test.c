#include <check.h>

extern TCase *linear_algebra_test();


int main(int argc, char **argv)
{
    Suite *s = suite_create("Core");
    suite_add_tcase(s, linear_algebra_test());

    SRunner *sr = srunner_create(s);

    srunner_run_all(sr, CK_ENV);
    int nf = srunner_ntests_failed(sr);
    srunner_free(sr);

    return nf == 0 ? 0 : 1;
}