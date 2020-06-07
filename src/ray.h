#ifndef RAY_H
#define RAY_H

struct Ray
{
    Vec3 origin;
    Vec3 direction;
    f32 time;
};

Point3 at(const Ray& ray, f32 t)
{
    return ray.origin + ray.direction * t;
}

Ray ray(Point3 origin, Vec3 direction, f32 time = 0.0f)
{
    Ray ray = {};
    ray.origin = origin;
    ray.direction = direction;
    ray.time = time;
    return ray;
}

#endif
