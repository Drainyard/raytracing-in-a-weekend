#ifndef RAY_H
#define RAY_H

struct Ray
{
    Vec3 origin;
    Vec3 direction;
};

Point3 at(const Ray& ray, f32 t)
{
    return ray.direction * t + ray.origin;
}

#endif
