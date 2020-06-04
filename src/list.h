#ifndef LIST_H
#define LIST_H

template<typename T>
struct List
{
    T *data;
    size_t count;

    size_t capacity;
};

template<typename T>
void maybe_grow(List<T>* list)
{
    if(list->count + 1 >= list->capacity)
    {
        if (list->capacity == 0)
        {
            list->capacity = 2;
        }
        else
        {
            list->capacity *= 2;
        }
        list->data = (T*)realloc(list->data, sizeof(T) * list->capacity);
    }
}

template<typename T>
T& add(List<T>* list, T value)
{
    maybe_grow(list);
    list->data[list->count++] = value;
    return list->data[list->count - 1];
}

template<typename T>
void clear(List<T>* list)
{
    if(list->capacity > 0)
    {
        list->capacity = 0;
        list->count = 0;
        free(list->data);
    }
}

#endif
