#include <check.h>

extern TCase *linear_algebra_test();
extern TCase *geometry_test();
extern TCase *kinematic_chain_test();


int main(int argc, char **argv)
{
    Suite *s = suite_create("Core");
    suite_add_tcase(s, linear_algebra_test());
    suite_add_tcase(s, geometry_test());
    suite_add_tcase(s, kinematic_chain_test());

    SRunner *sr = srunner_create(s);

    srunner_run_all(sr, CK_ENV);
    int nf = srunner_ntests_failed(sr);
    srunner_free(sr);

    return nf == 0 ? 0 : 1;
}