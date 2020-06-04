#ifndef MATH_UTIL_H
#define MATH_UTIL_H

#include <limits>

const f32 infinity = std::numeric_limits<f32>::infinity();
const f32 pi = 3.1415926535897932385;

#define RADS_IN_DEGREES (pi / 180.0f)

#define MIN(A,B) ((A) < (B) ? (A) : (B))
#define MAX(A,B) ((A) > (B) ? (A) : (B))

inline f32 degrees_to_radians(f32 degrees)
{
    return degrees * RADS_IN_DEGREES;
}

inline f32 fsqrt(f32 value)
{
    return (f32)sqrt(value);
}

inline f32 fcos(f32 value)
{
    return (f32)cos(value);
}

inline f32 fsin(f32 value)
{
    return (f32)sin(value);
}

inline f32 ftan(f32 value)
{
    return (f32)tan(value);
}

inline f32 fpow(f32 value, i32 exponent)
{
    return (f32)pow(value, exponent);
}

inline f32 random_float()
{
    return rand() / (RAND_MAX + 1.0f);
}

inline f32 random_float(f32 min, f32 max)
{
    return min + (max - min) * random_float();
}

inline f32 clamp(f32 x, f32 min, f32 max)
{
    if(x < min) return min;
    if(x > max) return max;
    return x;
}

#endif
