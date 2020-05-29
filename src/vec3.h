#ifndef VEC3_H
#define VEC3_H

struct Vec3
{
    union
    {
        struct
        {
            f32 x;
            f32 y;
            f32 z;
        };
        f32 e[3];
    };

    Vec3 operator-() const
    {
        return {-x, -y, -z};
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
    return {x, y, z};
}

Vec3 color(f32 x, f32 y, f32 z)
{
    return {x, y, z};
}

f32 length_squared(Vec3 v)
{
    return v[0] * v[0] + v[1] * v[1] + v[2] * v[2];
}

f32 length(Vec3 v)
{
    return (f32)sqrt(length_squared(v));
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

inline Vec3 unit_vector(Vec3 v)
{
    return v / length(v);
}

inline void print(const Vec3& v)
{
    printf("%f %f %f", v.x, v.y, v.z);
}



using Point3 = Vec3;
using Color = Vec3;

#endif
