#include <dyn2b/functions/mechanics.h>
#include <dyn2b/functions/linear_algebra.h>

#include <stdio.h>
#include <assert.h>


void mc_rbi_log(
        const struct mc_rbi *m)
{
    assert(m);

    printf("RigidBodyInertiaCoord(\n");
    printf("  I=                     h=       m=\n");
    for (int i = 0; i < 3; i++) {
        printf("  [%5.2f, %5.2f, %5.2f]",
                m->second_moment_of_mass.row[i].x,
                m->second_moment_of_mass.row[i].y,
                m->second_moment_of_mass.row[i].z);
        printf("  [%5.2f]", m->first_moment_of_mass.data[i]);
        if (i == 0) printf("  [%5.2f]", m->zeroth_moment_of_mass);
        if (i != 2) printf(",\n");
    }
    printf(")\n");
}


void ma_rbi_log(
        const struct ma_rbi *m)
{
    assert(m);
    assert(m->body);
    assert(m->point);
    assert(m->frame);

    printf("RigidBodyInertiaADT(point=%s|%s, frame={%s})\n",
        m->point->name,
        m->body->name,
        m->frame->name);
}