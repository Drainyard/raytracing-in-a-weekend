#ifndef ONB_H
#define ONB_H

union ONB
{
    Vec3 axis[3];
    struct
    {
        Vec3 u;
        Vec3 v;
        Vec3 w;
    };
    
    inline Vec3 operator[] (i32 i)
    {
        return axis[i];
    }
};

Vec3 local(ONB& onb, f32 a, f32 b, f32 c)
{
    return a * onb.u + b * onb.v + c * onb.w;
}

Vec3 local(ONB& onb, const Vec3& v)
{
    return v.x * onb.u + v.y * onb.v + v.z * onb.w;
}

ONB build_from_w(const Vec3& n)
{
    ONB onb = {};
    onb.axis[2] = unit_vector(n);
    Vec3 a = (fabs(onb.w.x) > 0.9f) ? vec3(0, 1, 0) : vec3(1, 0, 0);
    onb.axis[1] = unit_vector(cross(onb.w, a));
    onb.axis[0] = cross(onb.w, onb.v);
    return onb;
}


#endif
