#ifndef HITTABLE_H
#define HITTABLE_H

#include <algorithm>

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
    HITTABLE_FLIP_FACE,
    HITTABLE_BOX,
    HITTABLE_TRANSLATE,
    HITTABLE_ROTATE_Y,
    HITTABLE_CONSTANT_MEDIUM
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
        struct
        {
            Hittable* hittable;
        } flip_face;
        struct
        {
            List<Hittable> hittables;
            Point3 box_min;
            Point3 box_max;
        } box;
        struct
        {
            Hittable* hittable;
            Vec3 offset;
        } translate;
        struct
        {
            Hittable* hittable;
            f32 sin_theta;
            f32 cos_theta;
            b32 has_box;
            AABB bbox;
        } rotate_y;
        struct
        {
            Hittable* boundary;
            f32 neg_inv_density;
        } constant_medium;
    };
};

b32 hit(const List<Hittable>* list, const Ray& ray, f32 t_min, f32 t_max, Hit_Record& record);
b32 hit(const List<Hittable>* objects, const Hittable& hittable, const Ray& r, f32 t_min, f32 t_max, Hit_Record& record);

Hittable* alloc_hittable(Hittable hittable)
{
    Hittable* result = (Hittable*)malloc(sizeof(Hittable));
    *result = hittable;
    return result;
}

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
    Point3 small = point3(MIN(box0->min.x, box1->min.x),
                          MIN(box0->min.y, box1->min.y),
                          MIN(box0->min.z, box1->min.z));

    Point3 big = point3(MAX(box0->max.x, box1->max.x),
                        MAX(box0->max.y, box1->max.y),
                        MAX(box0->max.z, box1->max.z));

    return aabb(small, big);
}

bool bounding_box(List<Hittable>* list, const Hittable* hittable, f32 t0, f32 t1, AABB& output_box)
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
        f32 k = hittable->xy_rect.k;
        f32 x0 = hittable->xy_rect.x0;
        f32 x1 = hittable->xy_rect.x1;
        f32 y0 = hittable->xy_rect.y0;
        f32 y1 = hittable->xy_rect.y1;
        output_box = aabb(point3(x0, y0, k - 0.0001f),
                          point3(x1, y1, k + 0.0001f));
        return true;
    }
    break;
    case HITTABLE_XZ_RECT:
    {
        f32 x0 = hittable->xz_rect.x0;
        f32 x1 = hittable->xz_rect.x1;
        f32 z0 = hittable->xz_rect.z0;
        f32 z1 = hittable->xz_rect.z1;
        f32 k = hittable->xz_rect.k;
        output_box = aabb(point3(x0, k - 0.0001f, z0),
                          point3(x1, k + 0.0001f, z1));
        return true;
    }
    break;
    case HITTABLE_YZ_RECT:
    {
        f32 y0 = hittable->yz_rect.y0;
        f32 y1 = hittable->yz_rect.z1;
        f32 z0 = hittable->yz_rect.z0;
        f32 z1 = hittable->yz_rect.z1;
        f32 k = hittable->yz_rect.k;
        output_box = aabb(point3(k - 0.0001f, y0, z0),
                          point3(k + 0.0001f, y1, z1));
        return true;
    }
    break;
    case HITTABLE_FLIP_FACE:
    {
        return bounding_box(list, hittable->flip_face.hittable, t0, t1, output_box);
    }
    break;
    case HITTABLE_TRANSLATE:
    {
        if(!bounding_box(list, hittable->translate.hittable, t0, t1, output_box))
        {
            return false;
        }
        output_box = aabb(
            output_box.min + hittable->translate.offset,
            output_box.max + hittable->translate.offset
                          );

        return true;
    }
    break;
    case HITTABLE_BOX:
    {
        output_box = aabb(hittable->box.box_min, hittable->box.box_max);
        return true;
    }
    break;
    case HITTABLE_CONSTANT_MEDIUM:
    {
        return bounding_box(list, hittable->constant_medium.boundary, t0, t1, output_box);
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
        if(!bounding_box(list, &hittable, t0, t1, temp_box)) return false;
        output_box = first_box ? temp_box : surrounding_box(&output_box, &temp_box);
        first_box = false;
    }
    return true;
}

f32 pdf_value(const List<Hittable>* list, const Hittable& hittable, const Point3& origin, const Vec3& v)
{
    switch(hittable.type)
    {
    case HITTABLE_XZ_RECT:
    {
        Hit_Record rec = {};
        if (!hit(list, hittable, ray(origin, v), 0.001f, infinity, rec))
        {
            return 0.0f;
        }

        f32 area = (hittable.xz_rect.x1 - hittable.xz_rect.x0) * (hittable.xz_rect.z1 - hittable.xz_rect.z0);
        f32 distance_squared = rec.t * rec.t * length_squared(v);
        f32 cosine = fabs(dot(v, rec.normal) / length(v));

        return distance_squared / (cosine * area);
    }
    }
    return 0.0f;
}

Vec3 random(const Hittable& hittable, const Vec3& origin)
{
    switch(hittable.type)
    {
    case HITTABLE_XZ_RECT:
    {
        f32 x0 = hittable.xz_rect.x0;
        f32 x1 = hittable.xz_rect.x1;
        f32 z0 = hittable.xz_rect.z0;
        f32 z1 = hittable.xz_rect.z1;
        Point3 random_point = point3(random_float(x0, x1), hittable.xz_rect.k, random_float(z0, z1));
        return random_point - origin;
    }
    }
    return vec3(1.0f, 0.0f, 0.0f);
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

Hittable flip_face(Hittable* to_flip)
{
    Hittable hittable = {};
    hittable.type = HITTABLE_FLIP_FACE;
    hittable.flip_face.hittable = to_flip;
    return hittable;
}

Hittable box(Point3 p0, Point3 p1, size_t material)
{
    Hittable box = {};
    box.type = HITTABLE_BOX;

    box.box.box_min = p0;
    box.box.box_max = p1;

    add(&box.box.hittables, xy_rect(p0.x, p1.x, p0.y, p1.y, p1.z, material));
    add(&box.box.hittables, xy_rect(p0.x, p1.x, p0.y, p1.y, p0.z, material));
        
    add(&box.box.hittables, xz_rect(p0.x, p1.x, p0.z, p1.z, p1.y, material));
    add(&box.box.hittables, xz_rect(p0.x, p1.x, p0.z, p1.z, p0.y, material));
        
    add(&box.box.hittables, yz_rect(p0.y, p1.y, p0.z, p1.z, p1.x, material));
    add(&box.box.hittables, yz_rect(p0.y, p1.y, p0.z, p1.z, p0.x, material));

    return box;
}

Hittable translate(Hittable* hittable, Vec3 offset)
{
    Hittable translate = {};
    translate.type = HITTABLE_TRANSLATE;
    translate.translate.hittable = hittable;
    translate.translate.offset = offset;
    return translate;
}


Hittable rotate_y(List<Hittable>* list, Hittable* hittable, f32 angle)
{
    Hittable rotate_y = {};
    rotate_y.type = HITTABLE_ROTATE_Y;
    rotate_y.rotate_y.hittable = hittable;
    f32 radians = RADS_IN_DEGREES * angle;
    rotate_y.rotate_y.sin_theta = fsin(radians);
    rotate_y.rotate_y.cos_theta = fcos(radians);
    rotate_y.rotate_y.has_box = bounding_box(list, hittable, 0.0f, 1.0f, rotate_y.rotate_y.bbox);

    Point3 min = point3(infinity, infinity, infinity);
    Point3 max = point3(-infinity, -infinity, -infinity);

    f32 cos_theta = rotate_y.rotate_y.cos_theta;
    f32 sin_theta = rotate_y.rotate_y.sin_theta;

    AABB bbox = rotate_y.rotate_y.bbox;

    for(i32 i = 0; i < 2; i++)
    {
        for(i32 j = 0; j < 2; j++)
        {
            for(i32 k = 0; k < 2; k++)
            {
                f32 x = i * bbox.max.x + (1 - i) * bbox.min.x;
                f32 y = j * bbox.max.y + (1 - j) * bbox.min.y;
                f32 z = k * bbox.max.z + (1 - k) * bbox.min.z;

                f32 new_x = cos_theta * x + sin_theta * z;
                f32 new_z = -sin_theta * x + cos_theta * z;

                Vec3 tester = vec3(new_x, y, new_z);

                for(i32 c = 0; c < 3; c++)
                {
                    min[c] = MIN(min[c], tester[c]);
                    max[c] = MAX(max[c], tester[c]);
                }
            }
        }
    }

    bbox = aabb(min, max);
    
    return rotate_y;
}

Hittable constant_medium(List<Material>* material_list, Hittable* boundary, f32 density, size_t texture_handle)
{
    Hittable hittable = {};
    hittable.type = HITTABLE_CONSTANT_MEDIUM;
    hittable.material_handle = add(material_list, isotropic(texture_handle));

    hittable.constant_medium.boundary = boundary;
    hittable.constant_medium.neg_inv_density = -1.0f / density;

    return hittable;
}

List<Hittable>* temp_list;
inline i32 box_compare(List<Hittable>* list, const Hittable* a, const Hittable* b, i32 axis)
{
    temp_list = list;
    AABB box_a;
    AABB box_b;

    if(!bounding_box(list, a, 0.0f, 0.0f, box_a) || !bounding_box(list, b, 0.0f, 0.0f, box_b))
        printf("No bounding box in BVH_Node constructor\n");

    f32 val_a = box_a.min.e[axis];
    f32 val_b = box_b.min.e[axis];
    if(val_a < val_b)
    {
        return -1;
    }
    else if(val_a > val_b)
    {
        return 1;
    }

    return 0;
}

i32 box_x_compare(const void* a, const void* b)
{
    return box_compare(temp_list, (Hittable*)a, (Hittable*)b, 1);
}

i32 box_y_compare(const void* a, const void* b)
{
    return box_compare(temp_list, (Hittable*)a, (Hittable*)b, 1);
}

i32 box_z_compare(const void* a, const void* b)
{
    return box_compare(temp_list, (Hittable*)a, (Hittable*)b, 2);
}

Hittable bvh_node(List<Hittable>* objects, size_t start, size_t end, f32 time0, f32 time1)
{
    Hittable hittable = {};
    hittable.type = HITTABLE_BVH_NODE;
    temp_list = objects;
    
    i32 axis = random_int(0, 2);
    auto comparator = (axis == 0) ? box_x_compare :
        (axis == 1) ? box_y_compare : box_z_compare;

    size_t object_span = end - start;

    if(object_span == 1)
    {
        hittable.bvh_node.left = &objects->data[start];
        hittable.bvh_node.right = &objects->data[start];
    }
    else if(object_span == 2)
    {
        if(comparator(&objects->data[start], &objects->data[start + 1]) == -1)
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
        /* std::sort(objects->data + start, objects->data + end, comparator); */
        qsort(objects->data + start, end - start, sizeof(Hittable), comparator);

        size_t mid = start + (object_span/2);
        size_t left = add(objects, bvh_node(objects, start, mid, time0, time1));
        size_t right = add(objects, bvh_node(objects, mid, end, time0, time1));
        hittable.bvh_node.left = &objects->data[left];
        hittable.bvh_node.right = &objects->data[right];
    }
    
    AABB box_left, box_right;

    if(!bounding_box(objects, hittable.bvh_node.left, time0, time1, box_left) ||
       !bounding_box(objects, hittable.bvh_node.right, time0, time1, box_right))
    {
        printf("No bounding box in BVH_Node constructor\n");
    }

    hittable.bvh_node.box = surrounding_box(&box_left, &box_right);
    return hittable;
}

Hittable bvh_node(List<Hittable>* objects, f32 time0, f32 time1)
{
    return bvh_node(objects, 0, objects->count, time0, time1);
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

b32 hit(const List<Hittable>* objects, const Hittable& hittable, const Ray& r, f32 t_min, f32 t_max, Hit_Record& record)
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

        b32 hit_left = hit(objects, *hittable.bvh_node.left, r, t_min, t_max, record);
        b32 hit_right = hit(objects, *hittable.bvh_node.right, r, t_min, hit_left ? record.t : t_max, record);

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

        Vec3 outward_normal = vec3(0.0f, 1.0f, 0.0f);
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

        Vec3 outward_normal = vec3(1.0f, 0.0f, 0.0f);
        set_face_normal(record, r, outward_normal);
        record.material_handle = hittable.material_handle;
        record.p = at(r, t);
        return true;
    }
    break;
    case HITTABLE_FLIP_FACE:
    {
        if(!hit(objects, *hittable.flip_face.hittable, r, t_min, t_max, record))
        {
            return false;
        }

        record.front_face = !record.front_face;
        return true;
    }
    break;
    case HITTABLE_BOX:
    {
        return hit(&hittable.box.hittables, r, t_min, t_max, record);
    }
    break;
    case HITTABLE_TRANSLATE:
    {
        Ray moved_r = ray(r.origin - hittable.translate.offset, r.direction, r.time);
        if(!hit(objects, *hittable.translate.hittable, moved_r, t_min, t_max, record))
        {
            return false;
        }

        record.p += hittable.translate.offset;
        set_face_normal(record, moved_r, record.normal);
        return true;
    }
    break;
    case HITTABLE_ROTATE_Y:
    {
        Point3 origin = r.origin;
        Vec3 direction = r.direction;

        f32 cos_theta = hittable.rotate_y.cos_theta;
        f32 sin_theta = hittable.rotate_y.sin_theta;

        origin[0] = cos_theta * r.origin[0] - sin_theta * r.origin[2];
        origin[2] = sin_theta * r.origin[0] + cos_theta * r.origin[2];

        direction[0] = cos_theta * r.direction[0] - sin_theta * r.direction[2];
        direction[2] = sin_theta * r.direction[0] + cos_theta * r.direction[2];

        Ray rotated_r = ray(origin, direction, r.time);

        if(!hit(objects, *hittable.rotate_y.hittable, rotated_r, t_min, t_max, record))
        {
            return false;
        }

        Point3 p = record.p;
        Vec3 normal = record.normal;

        p[0] = cos_theta * record.p[0] + sin_theta * record.p[2];
        p[2] = -sin_theta * record.p[0] + cos_theta * record.p[2];

        normal[0] = cos_theta * record.normal[0] + sin_theta * record.normal[2];
        normal[2] = -sin_theta * record.normal[0] + cos_theta * record.normal[2];

        record.p = p;
        set_face_normal(record, rotated_r, normal);

        return true;
    }
    break;
    case HITTABLE_CONSTANT_MEDIUM:
    {
        const b32 enable_debug = false;
        const b32 debugging = enable_debug && random_float() < 0.00001;

        Hittable* boundary = hittable.constant_medium.boundary;        

        Hit_Record rec1 = {};
        Hit_Record rec2 = {};

        if(!hit(objects, *boundary, r, -infinity, infinity, rec1))
        {
            return false;
        }

        if(!hit(objects, *boundary, r, rec1.t + 0.0001f, infinity, rec2))
        {
            return false;
        }

        if(debugging)
        {
            printf("\nt0=%f, t1=%f", rec1.t, rec2.t);
        }

        if(rec1.t < t_min) rec1.t = t_min;
        if(rec2.t > t_max) rec2.t = t_max;

        if(rec1.t >= rec2.t)
        {
            return false;
        }

        if(rec1.t < 0)
        {
            rec1.t = 0.0f;
        }
        
        f32 ray_length = length(r.direction);
        f32 distance_inside_boundary = (rec2.t - rec1.t) * ray_length;
        f32 hit_distance = hittable.constant_medium.neg_inv_density * (f32)log(random_float());

        if(hit_distance > distance_inside_boundary)
        {
            return false;
        }

        record.t = rec1.t + hit_distance / ray_length;
        record.p = at(r, record.t);

        if(debugging)
        {
            printf("hit_distance = %f \n rec.t = %f\n rec.p = %f %f %f", hit_distance, record.t, record.p.x, record.p.y, record.p.z); 
        }

        record.normal = vec3(1, 0, 0);
        record.front_face = true;
        record.material_handle = hittable.material_handle;

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
        if(hit(list, hittable, ray, t_min, closest_so_far, temp_rec))
        {
            hit_anything = true;
            closest_so_far = temp_rec.t;
            record = temp_rec;
        }
    }
    return hit_anything;
}

#endif
