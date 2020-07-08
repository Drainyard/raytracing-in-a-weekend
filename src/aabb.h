#ifndef AABB_H
#define AABB_H

struct AABB
{
    Point3 min;
    Point3 max;
};

b32 hit(const AABB* aabb, const Ray& r, f32 tmin, f32 tmax)
{
    for(i32 a = 0; a < 3; a++)
    {
        f32 inv_d = 1.0f / r.direction[a];
        f32 t0 = (aabb->min[a] - r.origin[a]) * inv_d;
        f32 t1 = (aabb->max[a] - r.origin[a]) * inv_d;

        if(inv_d < 0.0f)
        {
            f32 t_temp = t0;
            t0 = t1;
            t1 = t_temp;
        }

        tmin = t0 > tmin ? t0 : tmin;
        tmax = t1 < tmax ? t1 : tmax;
        if(tmax <= tmin)
            return false;
    }
    return true;
}

AABB aabb(Point3 min, Point3 max)
{
    AABB result = {};
    result.min = min;
    result.max = max;
    return result;
}



#endif
