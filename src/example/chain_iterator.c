#include <dyn2b/example/chain_iterator.h>

#include <assert.h>


void kcc_iterator_reset(
        struct kcc_iterator_nbx *b)
{
    assert(b);
    assert(b->current_index);
    assert(b->next_index);
    assert(b->joint);
    assert(b->joint_type);
    assert(b->x_link);
    assert(b->has_next);

    int index = 0;

    *b->current_index = index;
    *b->next_index = index + 1;
    *b->joint = b->chain->segment[index].joint;
    *b->joint_type = b->chain->segment[index].joint.type;
    *b->x_link = b->chain->segment[index].joint_attachment;
    *b->has_next = (index < b->chain->number_of_segments);
}


void kcc_iterator_next(
        struct kcc_iterator_nbx *b)
{
    assert(b);
    assert(b->current_index);
    assert(b->next_index);
    assert(b->joint);
    assert(b->joint_type);
    assert(b->x_link);
    assert(b->has_next);

    int index = *b->current_index + 1;

    *b->current_index = index;
    *b->next_index = index + 1;
    *b->joint = b->chain->segment[index].joint;
    *b->joint_type = b->chain->segment[index].joint.type;
    *b->x_link = b->chain->segment[index].joint_attachment;
    *b->has_next = (index < b->chain->number_of_segments);
}
