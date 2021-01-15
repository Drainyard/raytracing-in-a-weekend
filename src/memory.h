#ifndef MEMORY_H
#define MEMORY_H

using umm = uintptr_t;

struct Block
{
    u8* base;
    u64 size;
    umm used;
    
    Block* prev;
    Block* next;
};

struct Allocator
{
    Block* current;
};

Allocator allocator;

inline umm get_alignment(umm alignment)
{
    umm alignment_offset = 0;

    umm ptr = (umm)allocator.current->base + allocator.current->used;
    umm mask = alignment - 1;
    if(ptr & mask)
    {
        alignment_offset = alignment - (ptr & mask);
    }
    return alignment_offset;
}

inline umm get_effective_size(size_t size_init)
{
    umm size = size_init;
    umm alignment_offset = get_alignment(8);

    size += alignment_offset;
    return size;
}

#define allocate_struct(type) (type*)allocate(sizeof(type))
void* allocate(size_t size_init)
{
    umm size = 0;

    if(allocator.current)
    {
        size = get_effective_size(size_init);
    }
    
    if(!allocator.current || (allocator.current->used + size) > allocator.current->size)
    {
        umm minimum_block_size = 1024 * 1024;
        u8* block = (u8*)malloc(minimum_block_size);
        Block* new_block = (Block*)malloc(sizeof(Block));
        new_block->base = block;
        new_block->prev = allocator.current;
        new_block->size = minimum_block_size;
        allocator.current = new_block;

        if(allocator.current->prev)
        {
            allocator.current->prev->next = new_block;
        }
    }

    umm alignment_offset = get_alignment(8);
    void* result = allocator.current->base + allocator.current->used + alignment_offset;
    allocator.current->used += size;

    return result;
}

void clear()
{
    Block* block = allocator.current;
    while(block)
    {
        Block* free_block = allocator.current;
        allocator.current = free_block->prev;
        free(free_block->base);
        free(free_block);
    }
}


#endif
