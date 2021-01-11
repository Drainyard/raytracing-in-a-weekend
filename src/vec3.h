#ifndef VEC3_H
#define VEC3_H

union Vec3
{
    struct
    {
        f32 x;
        f32 y;
        f32 z;
    };
    f32 e[3];

    Vec3 operator-() const
    {
        Vec3 result = {};
        result.x = -x;
        result.y = -y;
        result.z = -z;
        return result;
    }

    f32 operator[](i32 i) const
    {
        return e[i];
    }

    f32& operator[](i32 i)
    {
        return e[i];
    }

    Vec3& operator+=(const Vec3 &v)
    {
        e[0] += v.e[0];
        e[1] += v.e[1];
        e[2] += v.e[2];
        return *this;
    }

    Vec3& operator*=(const f32 t)
    {
        e[0] *= t;
        e[1] *= t;
        e[2] *= t;
        return *this;
    }

    Vec3& operator/=(const f32 t)
    {
        return *this *= 1/t;
    }
};

Vec3 vec3(f32 x, f32 y, f32 z)
{
    Vec3 result = {};
    result.e[0] = x;
    result.e[1] = y;
    result.e[2] = z;
    return result;
}

Vec3 point3(f32 x, f32 y, f32 z)
{
    return vec3(x, y, z);
}

Vec3 color(f32 x, f32 y, f32 z)
{
    return vec3(x, y, z);
}

f32 length_squared(Vec3 v)
{
    return v.e[0] * v.e[0] + v.e[1] * v.e[1] + v.e[2] * v.e[2];
}

f32 length(Vec3 v)
{
    return fsqrt(length_squared(v));
}

inline Vec3 operator+(const Vec3 &u, const Vec3 &v)
{
    return {u.e[0] + v.e[0], u.e[1] + v.e[1], u.e[2] + v.e[2]};
}

inline Vec3 operator-(const Vec3 &u, const Vec3 &v)
{
    return {u.e[0] - v.e[0], u.e[1] - v.e[1], u.e[2] - v.e[2]};
}

inline Vec3 operator*(const Vec3 &u, const Vec3 &v)
{
    return {u.e[0] * v.e[0], u.e[1] * v.e[1], u.e[2] * v.e[2]};
}

inline Vec3 operator*(f32 t, const Vec3 &v)
{
    return {t*v.e[0], t*v.e[1], t*v.e[2]};
}

inline Vec3 operator*(const Vec3 &v, f32 t)
{
    return t * v;
}

inline Vec3 operator/(Vec3 v, f32 t)
{
    return (1/t) * v;
}

inline f32 dot(const Vec3 &u, const Vec3 &v)
{
    return u.e[0] * v.e[0]
         + u.e[1] * v.e[1]
         + u.e[2] * v.e[2];
}

inline Vec3 cross(const Vec3 &u, const Vec3 &v)
{
    return {u.e[1] * v.e[2] - u.e[2] * v.e[1],
                u.e[2] * v.e[0] - u.e[0] * v.e[2],
            u.e[0] * v.e[1] - u.e[1] * v.e[0]};
}

inline bool near_zero(Vec3 v)
{
    const f32 s = 1e-8;
    return (fabs(v.e[0]) < s) && (fabs(v.e[1]) < s) && (fabs(v.e[2]) < s);
}

inline Vec3 unit_vector(Vec3 v)
{
    return v / length(v);
}

inline void print(const Vec3& v)
{
    printf("%f %f %f", v.x, v.y, v.z);
}

inline Vec3 random_vec3()
{
    Vec3 result = {};
    result.x = random_float();
    result.y = random_float();
    result.z = random_float();
    return result;
}

inline Vec3 random_vec3(f32 min, f32 max)
{
    Vec3 result = {};
    result.x = random_float(min, max);
    result.y = random_float(min, max);
    result.z = random_float(min, max);
    return result;
}

inline Vec3 random_in_unit_sphere()
{
    while(true)
    {
        Vec3 p = random_vec3(-1.0f, 1.0f);
        if(length_squared(p) >= 1.0f) continue;
        return p;
    }
}

inline Vec3 random_unit_vector()
{
    f32 a = random_float(0, 2.0f*pi);
    f32 z = random_float(-1.0f, 1.0f);
    f32 r = fsqrt(1 - z * z);
    Vec3 result = {};
    result.x = r * fcos(a);
    result.y = r * fsin(a);
    result.z = z;
    return result;
}

inline Vec3 random_in_hemisphere(const Vec3& normal)
{
    Vec3 result = random_in_unit_sphere();
    if(dot(result, normal) > 0.0f)
    {
        return result;
    }
    return -result;
}

inline Vec3 random_cosine_direction()
{
    f32 r1 = random_float();
    f32 r2 = random_float();
    f32 z = sqrt(1 - r2);

    f32 phi = 2 * pi * r1;
    f32 x = cos(phi) * sqrt(r2);
    f32 y = sin(phi) * sqrt(r2);

    return vec3(x, y, z);
}

inline Vec3 random_in_unit_disk()
{
    while(true)
    {
        Vec3 p = vec3(random_float(-1.0f, 1.0f), random_float(-1.0f, 1.0f), 0.0f);
        if(length_squared(p) >= 1.0f) continue;
        return p;
    }
}

inline Vec3 reflect(const Vec3& v, const Vec3& n)
{
    return v - 2 * dot(v, n) * n;
}

inline Vec3 refract(const Vec3& uv, const Vec3& n, f32 etai_over_etat)
{
    f32 cos_theta = dot(-uv, n);
    Vec3 r_out_parallel = etai_over_etat * (uv + cos_theta * n);
    Vec3 r_out_perp = -fsqrt(1.0f - length_squared(r_out_parallel)) * n;
    return r_out_parallel + r_out_perp;
}

using Point3 = Vec3;
using Color = Vec3;

#endif
