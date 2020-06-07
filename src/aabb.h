#ifndef AABB_H
#define AABB_H

struct AABB
{
    Point3 _min;
    Point3 _max;
};

b32 hit(const AABB* aabb, const Ray& r, f32 t_min, f32 t_max)
{
    for(i32 a = 0; a < 3; a++)
    {
        f32 inv_d = 1.0f / r.direction[a];
        f32 t0 = (aabb->_min[a] - r.origin[a]) * inv_d;
        f32 t1 = (aabb->_max[a] - r.origin[a]) * inv_d;

        if(inv_d < 0.0f)
        {
            f32 t_temp = t0;
            t0 = t1;
            t1 = t_temp;
        }

        t_min = t0 > t_min ? t0 : t_min;
        t_max = t1 < t_max ? t1 : t_max;
        if(t_max <= t_min)
            return false;
    }
    return true;
}

AABB aabb(Point3 min, Point3 max)
{
    AABB result = {};
    result._min = min;
    result._max = max;
    return result;
}



#endif
