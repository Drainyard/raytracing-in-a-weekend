#ifndef HITTABLE_H
#define HITTABLE_H

struct Hit_Record
{
    Point3 p;
    Vec3 normal;
    f32 t;
    b32 front_face;

    struct Material* material;
};

enum Hittable_Type
{
    HITTABLE_SPHERE
};

struct Hittable
{
    Hittable_Type type;
    struct Material* material;
    
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

Hittable sphere(Point3 center, f32 radius, Material* material)
{
    Hittable hittable = {};
    hittable.type = HITTABLE_SPHERE;
    hittable.sphere.center = center;
    hittable.sphere.radius = radius;
    hittable.material = material;
    return hittable;
}

inline void set_face_normal(Hit_Record& record, const Ray& r, const Vec3& outward_normal)
{
    record.front_face = dot(r.direction, outward_normal) < 0;
    record.normal = record.front_face ? outward_normal : -outward_normal;
}

b32 hit(const Hittable& hittable, const Ray& r, f32 t_min, f32 t_max, Hit_Record& record)
{
    switch(hittable.type)
    {
    case HITTABLE_SPHERE:
    {
        f32 radius = hittable.sphere.radius;
        Vec3 center = hittable.sphere.center;
        
        Vec3 oc = r.origin - center;
        f32 a = length_squared(r.direction);
        f32 half_b = dot(oc, r.direction);
        f32 c = length_squared(oc) - radius * radius;
        f32 discriminant = half_b * half_b - a * c;
        if(discriminant > 0)
        {
            f32 root = fsqrt(discriminant);
            f32 temp = (-half_b - root) / a;
            if(temp < t_max && temp > t_min)
            {
                record.t = temp;
                record.p = at(r, record.t);
                Vec3 outward_normal = (record.p - center) / radius;
                set_face_normal(record, r, outward_normal);
                record.material = hittable.material;
                return true;
            }
            temp = (-half_b + root) / a;
            if(temp < t_max && temp > t_min)
            {
                record.t = temp;
                record.p = at(r, record.t);
                Vec3 outward_normal = (record.p - center) / radius;
                set_face_normal(record, r, outward_normal);
                record.material = hittable.material;
                return true;
            }
        }
    }
    break;
    default:
    return false;
    }
    return false;
}


b32 hit(const Hittable_List* list, const Ray& ray, f32 t_min, f32 t_max, Hit_Record& record)
{
    Hit_Record temp_rec = {};
    b32 hit_anything = false;
    f32 closest_so_far = t_max;

    for(size_t i = 0; i < list->count; i++)
    {
        Hittable& hittable = list->hittables[i];
        if(hit(hittable, ray, t_min, closest_so_far, temp_rec))
        {
            hit_anything = true;
            closest_so_far = temp_rec.t;
            record = temp_rec;
        }
    }
    return hit_anything;
}

#endif
