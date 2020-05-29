#ifndef MATH_UTIL_H
#define MATH_UTIL_H

#include <limits>

const f32 infinity = std::numeric_limits<f32>::infinity();
const f32 pi = 3.1415926535897932385;

#define RADS_IN_DEGREES (pi / 180.0f)

inline f32 degrees_to_radians(f32 degrees)
{
    return degrees * RADS_IN_DEGREES;
}

inline f32 fsqrt(f32 value)
{
    return (f32)sqrt(value);
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
