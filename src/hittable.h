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
    HITTABLE_CONSTANT_MEDIUM,
    HITTABLE_LIST
};

struct Hittable
{
    Hittable_Type type;
    size_t material_handle;
};

struct Sphere : Hittable
{
    Point3 center;
    f32 radius;
};

struct MovingSphere : Hittable
{
    Point3 center0;
    Point3 center1;
    f32 radius;
    f32 time0;
    f32 time1;
};

struct BVHNode : Hittable
{
    Hittable* left;
    Hittable* right;
    AABB box;
};

struct XYRect : Hittable
{
    f32 x0;
    f32 x1;
    f32 y0;
    f32 y1;
    f32 k;
};

struct XZRect : Hittable
{
    f32 x0;
    f32 x1;
    f32 z0;
    f32 z1;
    f32 k;
};

struct YZRect : Hittable
{
    f32 y0;
    f32 y1;
    f32 z0;
    f32 z1;
    f32 k;
};

struct FlipFace : Hittable
{
    Hittable* ptr;
};

struct HittableList : Hittable
{
    List<Hittable*> list;
};

size_t add(HittableList* list, Hittable* value)
{
    return add(&list->list, value);
}

struct Box : Hittable
{
    Point3 box_min;
    Point3 box_max;
    HittableList sides;
};

struct Translate : Hittable
{
    Hittable* ptr;
    Vec3 offset;
};

struct RotateY : Hittable
{
    Hittable* ptr;
    f32 sin_theta;
    f32 cos_theta;
    b32 has_box;
    AABB bbox;
};

struct ConstantMedium : Hittable
{
    Hittable* boundary;
    f32 neg_inv_density;
};

#define AS(type, value) ((type*)value)


b32 hit(const HittableList& list, const Ray& ray, f32 t_min, f32 t_max, Hit_Record& record);
b32 hit(const HittableList& objects, const Hittable& hittable, const Ray& r, f32 t_min, f32 t_max, Hit_Record& record);

Hittable* alloc_hittable(Hittable hittable)
{
    Hittable* result = (Hittable*)malloc(sizeof(Hittable));
    *result = hittable;
    return result;
}

Point3 center(const MovingSphere* hittable, f32 time)
{
    Point3 center0 = hittable->center0;
    Point3 center1 = hittable->center1;
    f32 time0 = hittable->time0;
    f32 time1 = hittable->time1;
    
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

bool bounding_box(HittableList& list, const Hittable* hittable, f32 t0, f32 t1, AABB& output_box)
{
    switch(hittable->type)
    {
    case HITTABLE_SPHERE:
    {
        Sphere* sphere = AS(Sphere, hittable);
        f32 radius = sphere->radius;
        output_box = aabb(
            sphere->center - vec3(radius, radius, radius),
            sphere->center + vec3(radius, radius, radius)
                          );
        return true;
    }
    break;
    case HITTABLE_MOVING_SPHERE:
    {
        MovingSphere* moving_sphere = AS(MovingSphere, hittable);
        f32 radius = moving_sphere->radius;

        AABB box0 = aabb(center(moving_sphere, t0) - vec3(radius, radius, radius),
                         center(moving_sphere, t0) + vec3(radius, radius, radius));

        AABB box1 = aabb(center(moving_sphere, t1) - vec3(radius, radius, radius),
                         center(moving_sphere, t1) + vec3(radius, radius, radius));
        output_box = surrounding_box(&box0, &box1);
        return true;
    }
    break;
    case HITTABLE_BVH_NODE:
    {
        BVHNode* bvh = AS(BVHNode, hittable);
        output_box = bvh->box;
        return true;
    }
    break;
    case HITTABLE_XY_RECT:
    {
        XYRect* rect = AS(XYRect, hittable);
        f32 k = rect->k;
        f32 x0 = rect->x0;
        f32 x1 = rect->x1;
        f32 y0 = rect->y0;
        f32 y1 = rect->y1;
        output_box = aabb(point3(x0, y0, k - 0.0001f),
                          point3(x1, y1, k + 0.0001f));
        return true;
    }
    break;
    case HITTABLE_XZ_RECT:
    {
        XZRect* rect = AS(XZRect, hittable);
        f32 x0 = rect->x0;
        f32 x1 = rect->x1;
        f32 z0 = rect->z0;
        f32 z1 = rect->z1;
        f32 k = rect->k;
        output_box = aabb(point3(x0, k - 0.0001f, z0),
                          point3(x1, k + 0.0001f, z1));
        return true;
    }
    break;
    case HITTABLE_YZ_RECT:
    {
        YZRect* rect = AS(YZRect, hittable);
        f32 y0 = rect->y0;
        f32 y1 = rect->z1;
        f32 z0 = rect->z0;
        f32 z1 = rect->z1;
        f32 k = rect->k;
        output_box = aabb(point3(k - 0.0001f, y0, z0),
                          point3(k + 0.0001f, y1, z1));
        return true;
    }
    break;
    case HITTABLE_FLIP_FACE:
    {
        FlipFace* flip_face = AS(FlipFace, hittable);
        return bounding_box(list, flip_face->ptr, t0, t1, output_box);
    }
    break;
    case HITTABLE_TRANSLATE:
    {
        Translate* translate = AS(Translate, hittable);
        if(!bounding_box(list, translate->ptr, t0, t1, output_box))
        {
            return false;
        }
        output_box = aabb(
            output_box.min + translate->offset,
            output_box.max + translate->offset
                          );

        return true;
    }
    break;
    case HITTABLE_BOX:
    {
        Box* box = AS(Box, hittable);
        output_box = aabb(box->box_min, box->box_max);
        return true;
    }
    break;
    case HITTABLE_CONSTANT_MEDIUM:
    {
        ConstantMedium* medium = AS(ConstantMedium, hittable);
        return bounding_box(list, medium->boundary, t0, t1, output_box);
    }
    break;
    }
    return false;
}

bool bounding_box(HittableList& list, f32 t0, f32 t1, AABB& output_box)
{
    if(empty(&list.list)) return false;    

    AABB temp_box = {};
    b32 first_box = true;

    for(i32 i = 0; i < list.list.count; i++)
    {
        Hittable* hittable = list.list.data[i];
        if(!bounding_box(list, hittable, t0, t1, temp_box)) return false;
        output_box = first_box ? temp_box : surrounding_box(&output_box, &temp_box);
        first_box = false;
    }
    return true;
}

f32 pdf_value(const HittableList& list, const Hittable& hittable, const Point3& origin, const Vec3& v)
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
        
        XZRect* xz_rect = AS(XZRect, &hittable);

        f32 area = (xz_rect->x1 - xz_rect->x0) * (xz_rect->z1 - xz_rect->z0);
        f32 distance_squared = rec.t * rec.t * length_squared(v);
        f32 cosine = fabs(dot(v, rec.normal) / length(v));

        return distance_squared / (cosine * area);
    }
    case HITTABLE_SPHERE:
    {
        Hit_Record rec = {};
        if (!hit(list, hittable, ray(origin, v), 0.001f, infinity, rec))
        {
            return 0.0f;
        }

        Sphere* sphere = AS(Sphere, &hittable);
        f32 radius = sphere->radius;
        Point3 center = sphere->center;
        f32 cos_theta_max = sqrt(1 - radius * radius / length_squared(center - origin));
        f32 solid_angle = 2 * pi * (1 - cos_theta_max);

        return 1.0f / solid_angle;
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
        XZRect* xz_rect = AS(XZRect, &hittable);
        f32 x0 = xz_rect->x0;
        f32 x1 = xz_rect->x1;
        f32 z0 = xz_rect->z0;
        f32 z1 = xz_rect->z1;
        Point3 random_point = point3(random_float(x0, x1), xz_rect->k, random_float(z0, z1));
        return random_point - origin;
    }
    case HITTABLE_SPHERE:
    {
        Sphere* sphere = AS(Sphere, &hittable);
        Point3 center = sphere->center;
        f32 radius = sphere->radius;
        Vec3 direction = center - origin;
        f32 distance_squared = length_squared(direction);
        ONB uvw = build_from_w(direction);
        return local(uvw, random_to_sphere(radius, distance_squared));
    }
    }
    return vec3(1.0f, 0.0f, 0.0f);
}

Sphere* sphere(Point3 center, f32 radius, size_t material)
{
    Sphere* sphere = allocate_struct(Sphere);
    sphere->type = HITTABLE_SPHERE;
    sphere->material_handle = material;
    sphere->center = center;
    sphere->radius = radius;
    
    return sphere;
}

MovingSphere* moving_sphere(Point3 c0, Point3 c1, f32 t0, f32 t1, f32 radius, size_t material)
{
    MovingSphere* sphere = allocate_struct(MovingSphere);
    sphere->type = HITTABLE_MOVING_SPHERE;
    sphere->center0 = c0;
    sphere->center1 = c1;
    sphere->time0 = t0;
    sphere->time1 = t1;
    sphere->radius = radius;
    sphere->material_handle = material;
    return sphere;
}

XYRect* xy_rect(f32 x0, f32 x1, f32 y0, f32 y1, f32 k, size_t material)
{
    XYRect* xy_rect = allocate_struct(XYRect);
    xy_rect->type = HITTABLE_XY_RECT;
    xy_rect->material_handle = material;
    xy_rect->x0 = x0;
    xy_rect->x1 = x1;
    xy_rect->y0 = y0;
    xy_rect->y1 = y1;
    xy_rect->k = k;
    return xy_rect;
}

XZRect* xz_rect(f32 x0, f32 x1, f32 z0, f32 z1, f32 k, size_t material)
{
    XZRect* xz_rect = allocate_struct(XZRect);
    xz_rect->type = HITTABLE_XZ_RECT;
    xz_rect->material_handle = material;
    xz_rect->x0 = x0;
    xz_rect->x1 = x1;
    xz_rect->z0 = z0;
    xz_rect->z1 = z1;
    xz_rect->k = k;
    return xz_rect;
}

YZRect* yz_rect(f32 y0, f32 y1, f32 z0, f32 z1, f32 k, size_t material)
{
    YZRect* yz_rect = allocate_struct(YZRect);
    yz_rect->type = HITTABLE_YZ_RECT;
    yz_rect->material_handle = material;
    yz_rect->y0 = y0;
    yz_rect->y1 = y1;
    yz_rect->z0 = z0;
    yz_rect->z1 = z1;
    yz_rect->k = k;
    return yz_rect;
}

FlipFace* flip_face(Hittable* to_flip)
{
    FlipFace* flip_face = allocate_struct(FlipFace);
    flip_face->type = HITTABLE_FLIP_FACE;
    flip_face->ptr = to_flip;
    return flip_face;
}

Box* box(Point3 p0, Point3 p1, size_t material)
{
    Box* box = allocate_struct(Box);
    box->type = HITTABLE_BOX;

    box->box_min = p0;
    box->box_max = p1;

    add(&box->sides, xy_rect(p0.x, p1.x, p0.y, p1.y, p1.z, material));
    add(&box->sides, xy_rect(p0.x, p1.x, p0.y, p1.y, p0.z, material));
        
    add(&box->sides, xz_rect(p0.x, p1.x, p0.z, p1.z, p1.y, material));
    add(&box->sides, xz_rect(p0.x, p1.x, p0.z, p1.z, p0.y, material));
        
    add(&box->sides, yz_rect(p0.y, p1.y, p0.z, p1.z, p1.x, material));
    add(&box->sides, yz_rect(p0.y, p1.y, p0.z, p1.z, p0.x, material));

    return box;
}

Translate* translate(Hittable* hittable, Vec3 offset)
{
    Translate* translate = allocate_struct(Translate);
    translate->type = HITTABLE_TRANSLATE;
    translate->ptr = hittable;
    translate->offset = offset;
    return translate;
}


RotateY* rotate_y(HittableList& list, Hittable* hittable, f32 angle)
{
    RotateY* rotate_y = allocate_struct(RotateY);
    rotate_y->type = HITTABLE_ROTATE_Y;
    rotate_y->ptr = hittable;
    f32 radians = RADS_IN_DEGREES * angle;
    rotate_y->sin_theta = fsin(radians);
    rotate_y->cos_theta = fcos(radians);
    rotate_y->has_box = bounding_box(list, hittable, 0.0f, 1.0f, rotate_y->bbox);

    Point3 min = point3(infinity, infinity, infinity);
    Point3 max = point3(-infinity, -infinity, -infinity);

    f32 cos_theta = rotate_y->cos_theta;
    f32 sin_theta = rotate_y->sin_theta;

    AABB bbox = rotate_y->bbox;

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

ConstantMedium* constant_medium(List<Material>* material_list, Hittable* boundary, f32 density, size_t texture_handle)
{
    ConstantMedium* medium = allocate_struct(ConstantMedium);
    medium->type = HITTABLE_CONSTANT_MEDIUM;
    medium->material_handle = add(material_list, isotropic(texture_handle));

    medium->boundary = boundary;
    medium->neg_inv_density = -1.0f / density;

    return medium;
}

HittableList hittable_list()
{
    HittableList list = {};
    list.type = HITTABLE_LIST;
    return list;
}

HittableList* hittable_list(Hittable& object)
{
    HittableList* list = allocate_struct(HittableList);
    list->type = HITTABLE_LIST;
    add(list, &object);
    return list;
}

HittableList* temp_list;
inline i32 box_compare(HittableList* list, const Hittable* a, const Hittable* b, i32 axis)
{
    temp_list = list;
    AABB box_a;
    AABB box_b;

    if(!bounding_box(*list, a, 0.0f, 0.0f, box_a) || !bounding_box(*list, b, 0.0f, 0.0f, box_b))
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

BVHNode* bvh_node(HittableList& objects, size_t start, size_t end, f32 time0, f32 time1)
{
    BVHNode* bvh = allocate_struct(BVHNode);
    bvh->type = HITTABLE_BVH_NODE;
    temp_list = &objects;
    
    i32 axis = random_int(0, 2);
    auto comparator = (axis == 0) ? box_x_compare :
        (axis == 1) ? box_y_compare : box_z_compare;

    size_t object_span = end - start;

    if(object_span == 1)
    {
        bvh->left = objects.list.data[start];
        bvh->right = objects.list.data[start];
    }
    else if(object_span == 2)
    {
        if(comparator(&objects.list.data[start], &objects.list.data[start + 1]) == -1)
        {
            bvh->left = objects.list.data[start];
            bvh->right = objects.list.data[start + 1];
        }
        else
        {
            bvh->left = objects.list.data[start + 1];
            bvh->right = objects.list.data[start];
        }
    }
    else
    {
        /* std::sort(objects->data + start, objects->data + end, comparator); */
        qsort(objects.list.data + start, end - start, sizeof(Hittable), comparator);

        size_t mid = start + (object_span/2);
        size_t left = add(&objects, bvh_node(objects, start, mid, time0, time1));
        size_t right = add(&objects, bvh_node(objects, mid, end, time0, time1));
        bvh->left = objects.list.data[left];
        bvh->right = objects.list.data[right];
    }
    
    AABB box_left, box_right;

    if(!bounding_box(objects, bvh->left, time0, time1, box_left) ||
       !bounding_box(objects, bvh->right, time0, time1, box_right))
    {
        printf("No bounding box in BVH_Node constructor\n");
    }

    bvh->box = surrounding_box(&box_left, &box_right);
    return bvh;
}

BVHNode* bvh_node(HittableList& objects, f32 time0, f32 time1)
{
    return bvh_node(objects, 0, objects.list.count, time0, time1);
}

inline void set_face_normal(Hit_Record& record, const Ray& r, const Vec3& outward_normal)
{
    record.front_face = dot(r.direction, outward_normal) < 0;
    record.normal = record.front_face ? outward_normal : -outward_normal;
}

void get_sphere_uv(const Vec3& p, f32& u, f32& v)
{
    f32 theta = acos(-p.y);
    f32 phi = atan2(-p.z, p.x) + pi;

    u = phi / (2 * pi);
    v = theta / pi;
}

b32 hit(const HittableList& objects, const Hittable& hittable, const Ray& r, f32 t_min, f32 t_max, Hit_Record& record)
{
    switch(hittable.type)
    {
    case HITTABLE_SPHERE:
    {
        Sphere* sphere = AS(Sphere, &hittable);
        f32 radius = sphere->radius;
        Vec3 center = sphere->center;
        
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
        MovingSphere* sphere = AS(MovingSphere, &hittable);
        f32 radius = sphere->radius;
        Vec3 oc = r.origin - center(sphere, r.time);

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
                Vec3 cent = center(sphere, r.time);
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
                Vec3 cent = center(sphere, r.time);
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
        BVHNode* bvh = AS(BVHNode, &hittable);
        if(!hit(&bvh->box, r, t_min, t_max))
            return false;

        b32 hit_left = hit(objects, *bvh->left, r, t_min, t_max, record);
        b32 hit_right = hit(objects, *bvh->right, r, t_min, hit_left ? record.t : t_max, record);

        return hit_left || hit_right;
    }
    break;
    case HITTABLE_XY_RECT:
    {
        XYRect* xy_rect = AS(XYRect, &hittable);
        f32 x0 = xy_rect->x0;
        f32 x1 = xy_rect->x1;
        f32 y0 = xy_rect->y0;
        f32 y1 = xy_rect->y1;
        f32 k = xy_rect->k;
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
        XZRect* xz_rect = AS(XZRect, &hittable);
        f32 x0 = xz_rect->x0;
        f32 x1 = xz_rect->x1;
        f32 z0 = xz_rect->z0;
        f32 z1 = xz_rect->z1;
        f32 k = xz_rect->k;
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
        YZRect* yz_rect = AS(YZRect, &hittable);
        f32 y0 = yz_rect->y0;
        f32 y1 = yz_rect->y1;
        f32 z0 = yz_rect->z0;
        f32 z1 = yz_rect->z1;
        f32 k =  yz_rect->k;
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
        FlipFace* flip_face = AS(FlipFace, &hittable);
        if(!hit(objects, *flip_face->ptr, r, t_min, t_max, record))
        {
            return false;
        }

        record.front_face = !record.front_face;
        return true;
    }
    break;
    case HITTABLE_BOX:
    {
        Box* box = AS(Box, &hittable);
        return hit(box->sides, r, t_min, t_max, record);
    }
    break;
    case HITTABLE_TRANSLATE:
    {
        Translate* translate = AS(Translate, &hittable);
        Ray moved_r = ray(r.origin - translate->offset, r.direction, r.time);
        if(!hit(objects, *translate->ptr, moved_r, t_min, t_max, record))
        {
            return false;
        }

        record.p += translate->offset;
        set_face_normal(record, moved_r, record.normal);
        return true;
    }
    break;
    case HITTABLE_ROTATE_Y:
    {
        RotateY* rotate_y = AS(RotateY, &hittable);
        Point3 origin = r.origin;
        Vec3 direction = r.direction;

        f32 cos_theta = rotate_y->cos_theta;
        f32 sin_theta = rotate_y->sin_theta;

        origin[0] = cos_theta * r.origin[0] - sin_theta * r.origin[2];
        origin[2] = sin_theta * r.origin[0] + cos_theta * r.origin[2];

        direction[0] = cos_theta * r.direction[0] - sin_theta * r.direction[2];
        direction[2] = sin_theta * r.direction[0] + cos_theta * r.direction[2];

        Ray rotated_r = ray(origin, direction, r.time);

        if(!hit(objects, *rotate_y->ptr, rotated_r, t_min, t_max, record))
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

        ConstantMedium* medium = AS(ConstantMedium, &hittable);

        Hittable* boundary = medium->boundary;        

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
        f32 hit_distance = medium->neg_inv_density * (f32)log(random_float());

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


b32 hit(const HittableList& list, const Ray& ray, f32 t_min, f32 t_max, Hit_Record& record)
{
    Hit_Record temp_rec = {};
    b32 hit_anything = false;
    f32 closest_so_far = t_max;

    for(size_t i = 0; i < list.list.count; i++)
    {
        Hittable* hittable = list.list.data[i];
        if(hit(list, *hittable, ray, t_min, closest_so_far, temp_rec))
        {
            hit_anything = true;
            closest_so_far = temp_rec.t;
            record = temp_rec;
        }
    }
    return hit_anything;
}

#endif
