#ifndef HITTABLE_H
#define HITTABLE_H

struct Hit_Record
{
    Point3 p;
    Vec3 normal;
    f32 t;
    b32 front_face;
    f32 u;
    f32 v;

    size_t material_handle;
};

enum Hittable_Type
{
    HITTABLE_SPHERE,
    HITTABLE_MOVING_SPHERE,
    HITTABLE_BVH_NODE,
    HITTABLE_XY_RECT,
    HITTABLE_XZ_RECT,
    HITTABLE_YZ_RECT,
};

struct Hittable
{
    Hittable_Type type;
    size_t material_handle;
    
    union
    {
        struct
        {
            Point3 center;
            f32 radius;
        } sphere;
        struct
        {
            Point3 center0;
            Point3 center1;
            f32 time0;
            f32 time1;
            f32 radius;
        } moving_sphere;
        struct
        {
            Hittable* left;
            Hittable* right;
            AABB box;
        } bvh_node;
        struct
        {
            f32 x0;
            f32 x1;
            f32 y0;
            f32 y1;
            f32 k;
        } xy_rect;
        struct
        {
            f32 x0;
            f32 x1;
            f32 z0;
            f32 z1;
            f32 k;
        } xz_rect;
        struct
        {
            f32 y0;
            f32 y1;
            f32 z0;
            f32 z1;
            f32 k;
        } yz_rect;
    };
};

Point3 center(const Hittable* hittable, f32 time)
{
    if(hittable->type != HITTABLE_MOVING_SPHERE)
    {
        return point3(0, 0, 0);
    }
    
    Point3 center0 = hittable->moving_sphere.center0;
    Point3 center1 = hittable->moving_sphere.center1;
    f32 time0 = hittable->moving_sphere.time0;
    f32 time1 = hittable->moving_sphere.time1;
    
    return center0 + ((time - time0) / (time1 - time0)) * (center1 - center0);
}

AABB surrounding_box(AABB *box0, AABB* box1)
{
    Point3 small = point3(MIN(box0->_min.x, box1->_min.x),
                          MIN(box0->_min.y, box1->_min.y),
                          MIN(box0->_min.z, box1->_min.z));

    Point3 big = point3(MAX(box0->_max.x, box1->_max.x),
                        MAX(box0->_max.y, box1->_max.y),
                        MAX(box0->_max.z, box1->_max.z));

    return aabb(small, big);
}

bool bounding_box(const Hittable* hittable, f32 t0, f32 t1, AABB& output_box)
{
    switch(hittable->type)
    {
    case HITTABLE_SPHERE:
    {
        f32 radius = hittable->sphere.radius;
        output_box = aabb(
            hittable->sphere.center - vec3(radius, radius, radius),
            hittable->sphere.center + vec3(radius, radius, radius)
                          );
        return true;
    }
    break;
    case HITTABLE_MOVING_SPHERE:
    {
        f32 radius = hittable->moving_sphere.radius;
        output_box = aabb(
            hittable->sphere.center - vec3(radius, radius, radius),
            hittable->sphere.center + vec3(radius, radius, radius)
                          );
        AABB box0 = aabb(center(hittable, t0) - vec3(radius, radius, radius),
                         center(hittable, t0) + vec3(radius, radius, radius));

        AABB box1 = aabb(center(hittable, t1) - vec3(radius, radius, radius),
                         center(hittable, t1) + vec3(radius, radius, radius));
        output_box = surrounding_box(&box0, &box1);
        return true;
    }
    break;
    case HITTABLE_BVH_NODE:
    {
        output_box = hittable->bvh_node.box;
        return true;
    }
    break;
    case HITTABLE_XY_RECT:
    {
        output_box = aabb(point3(hittable->xy_rect.x0, hittable->xy_rect.y0, hittable->xy_rect.k - 0.0001f),
                          point3(hittable->xy_rect.x1, hittable->xy_rect.y1, hittable->xy_rect.k + 0.0001f));
        return true;
    }
    break;
    case HITTABLE_XZ_RECT:
    {
        output_box = aabb(point3(hittable->xz_rect.x0, hittable->xz_rect.k - 0.0001f, hittable->xz_rect.z0),
                          point3(hittable->xz_rect.x1, hittable->xz_rect.k + 0.0001f, hittable->xz_rect.z1));
        return true;
    }
    break;
    case HITTABLE_YZ_RECT:
    {
        output_box = aabb(point3(hittable->yz_rect.k - 0.0001f, hittable->yz_rect.y0, hittable->yz_rect.z0),
                          point3(hittable->yz_rect.k + 0.0001f, hittable->yz_rect.y1, hittable->yz_rect.z1));
        return true;
    }
    break;
    }
    return false;
}

bool bounding_box(List<Hittable>* list, f32 t0, f32 t1, AABB& output_box)
{
    if(empty(list)) return false;

    AABB temp_box = {};
    b32 first_box = true;

    for(i32 i = 0; i < list->count; i++)
    {
        Hittable& hittable = list->data[i];
        if(!bounding_box(&hittable, t0, t1, temp_box)) return false;
        output_box = first_box ? temp_box : surrounding_box(&output_box, &temp_box);
        first_box = false;
    }
    return true;
}


Hittable sphere(Point3 center, f32 radius, size_t material)
{
    Hittable hittable = {};
    hittable.type = HITTABLE_SPHERE;
    hittable.sphere.center = center;
    hittable.sphere.radius = radius;
    hittable.material_handle = material;
    return hittable;
}

Hittable moving_sphere(Point3 c0, Point3 c1, f32 t0, f32 t1, f32 radius, size_t material)
{
    Hittable hittable = {};
    hittable.type = HITTABLE_MOVING_SPHERE;
    hittable.moving_sphere.center0 = c0;
    hittable.moving_sphere.center1 = c1;
    hittable.moving_sphere.time0 = t0;
    hittable.moving_sphere.time1 = t1;
    hittable.moving_sphere.radius = radius;
    hittable.material_handle = material;
    return hittable;
}

Hittable xy_rect(f32 x0, f32 x1, f32 y0, f32 y1, f32 k, size_t material)
{
    Hittable hittable = {};
    hittable.type = HITTABLE_XY_RECT;
    hittable.material_handle = material;
    hittable.xy_rect.x0 = x0;
    hittable.xy_rect.x1 = x1;
    hittable.xy_rect.y0 = y0;
    hittable.xy_rect.y1 = y1;
    hittable.xy_rect.k = k;
    return hittable;
}

Hittable xz_rect(f32 x0, f32 x1, f32 z0, f32 z1, f32 k, size_t material)
{
    Hittable hittable = {};
    hittable.type = HITTABLE_XZ_RECT;
    hittable.material_handle = material;
    hittable.xz_rect.x0 = x0;
    hittable.xz_rect.x1 = x1;
    hittable.xz_rect.z0 = z0;
    hittable.xz_rect.z1 = z1;
    hittable.xz_rect.k = k;
    return hittable;
}

Hittable yz_rect(f32 y0, f32 y1, f32 z0, f32 z1, f32 k, size_t material)
{
    Hittable hittable = {};
    hittable.type = HITTABLE_YZ_RECT;
    hittable.material_handle = material;
    hittable.yz_rect.y0 = y0;
    hittable.yz_rect.y1 = y1;
    hittable.yz_rect.z0 = z0;
    hittable.yz_rect.z1 = z1;
    hittable.yz_rect.k = k;
    return hittable;
}

inline i32 box_compare(const Hittable* a, const Hittable* b, i32 axis)
{
    AABB box_a;
    AABB box_b;

    if(!bounding_box(a, 0.0f, 0.0f, box_a) || !bounding_box(b, 0.0f, 0.0f, box_b))
        printf("No bounding box in BVH_Node constructor\n");

    return box_a._min.e[axis] < box_b._min.e[axis];
}

i32 box_x_compare(const void* a, const void* b)
{
    return box_compare((Hittable*)a, (Hittable*)b, 1);
}

i32 box_y_compare(const void* a, const void* b)
{
    return box_compare((Hittable*)a, (Hittable*)b, 1);
}

i32 box_z_compare(const void* a, const void* b)
{
    return box_compare((Hittable*)a, (Hittable*)b, 2);
}

Hittable bvh_node(List<Hittable>* objects, size_t start, size_t end, f32 time0, f32 time1)
{
    Hittable hittable = {};
    hittable.type = HITTABLE_BVH_NODE;
    
    i32 axis = random_int(0, 2);
    i32 (*comparator)(const void*, const void*) = (axis == 0) ? box_x_compare :
        (axis == 1) ? box_y_compare : box_z_compare;

    size_t object_span = end - start;

    if(object_span == 1)
    {
        hittable.bvh_node.left = &objects->data[start];
    }
    else if(object_span == 2)
    {
        if(comparator(&objects->data[start], &objects->data[start + 1]))
        {
            hittable.bvh_node.left = &objects->data[start];
            hittable.bvh_node.right = &objects->data[start + 1];
        }
        else
        {
            hittable.bvh_node.left = &objects->data[start + 1];
            hittable.bvh_node.right = &objects->data[start];
        }
    }
    else
    {
        qsort(objects->data, objects->count, sizeof(Hittable), comparator);

        size_t mid = start + (object_span/2);
        Hittable left = bvh_node(objects, start, mid, time0, time1);
        Hittable right = bvh_node(objects, mid, end, time0, time1);
        hittable.bvh_node.left = (Hittable*)malloc(sizeof(Hittable));
        hittable.bvh_node.right = (Hittable*)malloc(sizeof(Hittable));
        *hittable.bvh_node.left = left;
        *hittable.bvh_node.right = right;
    }

    AABB box_left, box_right;

    if(!bounding_box(hittable.bvh_node.left, time0, time1, box_left) ||
       !bounding_box(hittable.bvh_node.right, time0, time1, box_right))
    {
        printf("No bounding box in BVH_Node constructor\n");
    }

    hittable.bvh_node.box = surrounding_box(&box_left, &box_right);
    return hittable;
}

inline void set_face_normal(Hit_Record& record, const Ray& r, const Vec3& outward_normal)
{
    record.front_face = dot(r.direction, outward_normal) < 0;
    record.normal = record.front_face ? outward_normal : -outward_normal;
}

void get_sphere_uv(const Vec3& p, f32& u, f32& v)
{
    f32 phi = (f32)atan2(p.z, p.x);
    f32 theta = (f32)asin(p.y);
    u = 1.0f - (phi + pi) / (2 * pi);
    v = (theta + pi / 2) / pi;
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
                get_sphere_uv((record.p - center) / radius, record.u, record.v);
                record.material_handle = hittable.material_handle;
                return true;
            }
            temp = (-half_b + root) / a;
            if(temp < t_max && temp > t_min)
            {
                record.t = temp;
                record.p = at(r, record.t);
                Vec3 outward_normal = (record.p - center) / radius;
                set_face_normal(record, r, outward_normal);
                get_sphere_uv((record.p - center) / radius, record.u, record.v);
                record.material_handle = hittable.material_handle;
                return true;
            }
        }
    }
    break;
    case HITTABLE_MOVING_SPHERE:
    {
        f32 radius = hittable.moving_sphere.radius;
        Vec3 oc = r.origin - center(&hittable, r.time);

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
                Vec3 cent = center(&hittable, r.time);
                Vec3 p_minus_c = record.p - cent;
                Vec3 outward_normal = (p_minus_c) / radius;
                set_face_normal(record, r, outward_normal);
                get_sphere_uv(p_minus_c / radius, record.u, record.v);
                record.material_handle = hittable.material_handle;
                return true;
            }
            temp = (-half_b + root) / a;
            if(temp < t_max && temp > t_min)
            {
                record.t = temp;
                record.p = at(r, record.t);
                Vec3 cent = center(&hittable, r.time);
                Vec3 p_minus_c = record.p - cent;
                Vec3 outward_normal = (p_minus_c) / radius;
                set_face_normal(record, r, outward_normal);
                get_sphere_uv(p_minus_c / radius, record.u, record.v);
                record.material_handle = hittable.material_handle;
                return true;
            }
        }
    }
    break;
    case HITTABLE_BVH_NODE:
    {
        if(!hit(&hittable.bvh_node.box, r, t_min, t_max))
            return false;

        b32 hit_left = hit(*hittable.bvh_node.left, r, t_min, t_max, record);
        b32 hit_right = hit(*hittable.bvh_node.right, r, t_min, hit_left ? record.t : t_max, record);

        return hit_left || hit_right;
    }
    break;
    case HITTABLE_XY_RECT:
    {
        f32 x0 = hittable.xy_rect.x0;
        f32 x1 = hittable.xy_rect.x1;
        f32 y0 = hittable.xy_rect.y0;
        f32 y1 = hittable.xy_rect.y1;
        f32 k = hittable.xy_rect.k;
        // Solve for t:
        f32 t = (k - r.origin.z) / r.direction.z;
        if(t < t_min || t > t_max)
        {
            return false;
        }

        f32 x = r.origin.x + r.direction.x * t;
        f32 y = r.origin.y + r.direction.y * t;

        if(x < x0 || x > x1 || y < y0 || y > y1)
        {
            return false;
        }

        record.u = (x - x0) / (x1 - x0);
        record.v = (y - y0) / (y1 - y0);
        record.t = t;

        Vec3 outward_normal = vec3(0.0f, 0.0f, 1.0f);
        set_face_normal(record, r, outward_normal);
        record.material_handle = hittable.material_handle;
        record.p = at(r, t);
        return true;
    }
    break;
    case HITTABLE_XZ_RECT:
    {
        f32 x0 = hittable.xz_rect.x0;
        f32 x1 = hittable.xz_rect.x1;
        f32 z0 = hittable.xz_rect.z0;
        f32 z1 = hittable.xz_rect.z1;
        f32 k = hittable.xz_rect.k;
        // Solve for t:
        f32 t = (k - r.origin.y) / r.direction.y;
        if(t < t_min || t > t_max)
        {
            return false;
        }

        f32 x = r.origin.x + r.direction.x * t;
        f32 z = r.origin.z + r.direction.z * t;

        if(x < x0 || x > x1 || z < z0 || z > z1)
        {
            return false;
        }

        record.u = (x - x0) / (x1 - x0);
        record.v = (z - z0) / (z1 - z0);
        record.t = t;

        Vec3 outward_normal = vec3(0.0f, 0.0f, 1.0f);
        set_face_normal(record, r, outward_normal);
        record.material_handle = hittable.material_handle;
        record.p = at(r, t);
        return true;
    }
    break;
    case HITTABLE_YZ_RECT:
    {
        f32 y0 = hittable.yz_rect.y0;
        f32 y1 = hittable.yz_rect.y1;
        f32 z0 = hittable.yz_rect.z0;
        f32 z1 = hittable.yz_rect.z1;
        f32 k = hittable.yz_rect.k;
        // Solve for t:
        f32 t = (k - r.origin.x) / r.direction.x;
        if(t < t_min || t > t_max)
        {
            return false;
        }

        f32 y = r.origin.y + r.direction.y * t;
        f32 z = r.origin.z + r.direction.z * t;

        if(y < y0 || y > y1 || z < z0 || z > z1)
        {
            return false;
        }

        record.u = (y - y0) / (y1 - y0);
        record.v = (z - z0) / (z1 - z0);
        record.t = t;

        Vec3 outward_normal = vec3(0.0f, 0.0f, 1.0f);
        set_face_normal(record, r, outward_normal);
        record.material_handle = hittable.material_handle;
        record.p = at(r, t);
        return true;
    }
    break;
    default:
    return false;
    }
    return false;
}


b32 hit(const List<Hittable>* list, const Ray& ray, f32 t_min, f32 t_max, Hit_Record& record)
{
    Hit_Record temp_rec = {};
    b32 hit_anything = false;
    f32 closest_so_far = t_max;

    for(size_t i = 0; i < list->count; i++)
    {
        Hittable& hittable = list->data[i];
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
