#include <dyn2b/functions/kinematic_chain_nbx.h>
#include <dyn2b/functions/kinematic_chain.h>

#include <assert.h>


void kcc_fpk_nbx(
        struct kcc_fpk_nbx *nbx)
{
    assert(nbx);
    assert(nbx->joint);

    kcc_joint[nbx->joint->type].fpk(nbx->joint, nbx->q, nbx->x);
}
