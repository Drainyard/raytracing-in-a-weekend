#ifndef HITTABLE_H
#define HITTABLE_H

struct Hit_Record
{
    Point3 p;
    Vec3 normal;
    f32 t;
    b32 front_face;
};

enum Hittable_Type
{
    HITTABLE_SPHERE
};

struct Hittable
{
    Hittable_Type type;

    union
    {
        struct
        {
            Point3 center;
            f32 radius;
        } sphere;
    };
};

struct Hittable_List
{
    Hittable *hittables;
    size_t count;

    size_t capacity;
};

void maybe_grow(Hittable_List* list)
{
    if(list->count + 1 == list->capacity)
    {
        list->capacity *= 2;
        list->hittables = (Hittable*)realloc(list->hittables, sizeof(Hittable) * list->capacity);
    }
}

void add(Hittable_List* list, Hittable hittable)
{
    maybe_grow(list);
    list->hittables[list->count++] = hittable;
}

void clear(Hittable_List* list)
{
    if(list->capacity > 0)
    {
        list->capacity = 0;
        list->count = 0;
        free(list->hittables);
    }
}

inline void set_face_normal(Hit_Record& record, const Ray& r, const Vec3& outward_normal)
{
    record.front_face = dot(r.direction, outward_normal) < 0;
    record.normal = record.front_face ? outward_normal : -outward_normal;
}

b32 hit(Hittable& hittable, const Ray& r, f32 t_min, f32 t_max, Hit_Record& record)
{
    switch(hittable.type)
    {
    case HITTABLE_SPHERE:
    {
        Vec3 oc = r.origin - hittable.sphere.center;
        f32 a = length_squared(r.direction);
        f32 half_b = dot(oc, r.direction);
        f32 c = length_squared(oc) - hittable.sphere.radius * hittable.sphere.radius;
        f32 discriminant = half_b * half_b - a * c;
        if(discriminant > 0)
        {
            f32 root = sqrt(discriminant);
            f32 temp = (-half_b - root) / a;
            if(temp < t_max && temp > t_min)
            {
                record.t = temp;
                record.p = at(r, record.t);
                Vec3 outward_normal = (record.p - hittable.sphere.center) / hittable.sphere.radius;
                set_face_normal(record, r, outward_normal);
                return true;
            }
        }
        temp = (-half_b + root) / a;
        if(temp < t_max && temp > t_min)
        {
            record.t = temp;
            record.p = at(r, record.t);
            Vec3 outward_normal = (record.p - hittable.sphere.center) / hittable.sphere.radius;
            set_face_normal(record, r, outward_normal);
            return true;
        }
    }
    break;
    default:
    return false;
    }
    return false;
}


b32 hit(Hittable_List* list, const Ray& ray, f32 t_min, f32 t_max, Hit_Record& record)
{
    Hit_Record temp_rec = {};
    b32 hit_anything = false;
    f32 closest_so_far = t_max;

    for(size_t i = 0; i < list->count; i++)
    {
        Hittable& hittable = list->hittables[i];
        if(hit(hittable, ray, t_min, t_max, record))
        {
            hit_anything = true;
            closest_so_far = temp_rec.t;
            rec = temp_rec;
        }
    }
    return hit_anything;
}

#endif