#ifndef EXAMPLE_CHAIN_ITERATOR_H
#define EXAMPLE_CHAIN_ITERATOR_H

#include <dyn2b/types/kinematic_chain.h>

#include <stdbool.h>

#ifdef __cplusplus
extern "C" {
#endif


struct kcc_iterator_nbx
{
    struct kcc_kinematic_chain *chain;
    int *current_index;
    int *next_index;
    struct kcc_joint *joint;
    enum joint_type *joint_type;
    struct gc_pose *x_link;
    bool *has_next;
};


void kcc_iterator_reset(
        struct kcc_iterator_nbx *b);

void kcc_iterator_next(
        struct kcc_iterator_nbx *b);

#ifdef __cplusplus
}
#endif

#endif
