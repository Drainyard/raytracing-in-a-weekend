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

#endif
