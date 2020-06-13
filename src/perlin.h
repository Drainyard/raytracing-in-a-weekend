#ifndef PERLIN_H
#define PERLIN_H

static const int POINT_COUNT = 256;
struct Perlin
{
    Vec3* rand_vec;
    i32* perm_x;
    i32* perm_y;
    i32* perm_z;
};

void permute(i32* p, i32 n)
{
    for(i32 i = n - 1; i > 0; i--)
    {
        i32 target = random_int(0, i);
        i32 tmp = p[i];
        p[i] = p[target];
        p[target] = tmp;
    }
}

i32* perlin_generate_perm()
{
    i32* p = (i32*)malloc(sizeof(i32) * POINT_COUNT);

    for(i32 i = 0; i < POINT_COUNT; i++)
    {
        p[i] = i;
    }

    permute(p, POINT_COUNT);
    
    return p;
}

Perlin perlin()
{
    Perlin perlin = {};
    perlin.rand_vec = (Vec3*)malloc(sizeof(Vec3) * POINT_COUNT);
    for(i32 i = 0; i < POINT_COUNT; i++)
    {
        perlin.rand_vec[i] = unit_vector(random_vec3(-1.0f, 1.0f));
    }

    perlin.perm_x = perlin_generate_perm();
    perlin.perm_y = perlin_generate_perm();
    perlin.perm_z = perlin_generate_perm();
    
    return perlin;                  
}

inline f32 perlin_interp(Vec3 c[2][2][2], f32 u, f32 v, f32 w)
{
    f32 uu = u * u * (3 - 2 * u);
    f32 vv = v * v * (3 - 2 * v);
    f32 ww = w * w * (3 - 2 * w);
    f32 accum = 0.0f;
    
    for(i32 i = 0; i < 2; i++)
    {
        for(i32 j = 0; j < 2; j++)
        {
            for(i32 k = 0; k < 2; k++)
            {
                Vec3 weight_v = vec3(u - i, v - j, w - k);
                accum += (i * uu + (1 - i) * (1 - uu))
                    * (j * vv + (1 - j) * (1 - vv))
                    * (k * ww + (1 - k) * (1 - ww))
                    * dot(c[i][j][k], weight_v);
            }
        }
    }
    return accum;
}

f32 noise(Perlin* perlin, const Point3& p)
{
    f32 u = p.x - ffloor(p.x);
    f32 v = p.y - ffloor(p.y);
    f32 w = p.z - ffloor(p.z);

    i32 i = (i32)ffloor(p.x);
    i32 j = (i32)ffloor(p.y);
    i32 k = (i32)ffloor(p.z);
    Vec3 c[2][2][2];

    for(i32 di = 0; di < 2; di++)
    {
        for(i32 dj = 0; dj < 2; dj++)
        {
            for(i32 dk = 0; dk < 2; dk++)
            {
                c[di][dj][dk] = perlin->rand_vec[
                    perlin->perm_x[(i + di) & 255] ^
                    perlin->perm_y[(j + dj) & 255] ^
                    perlin->perm_z[(k + dk) & 255]
                                                   ];
            }
        }
    }

    return perlin_interp(c, u, v, w);
}

f32 turb(Perlin* perlin, const Point3& p, i32 depth = 7)
{
    f32 accum = 0.0f;
    Point3 temp_p = p;
    f32 weight = 1.0f;

    for(i32 i = 0; i < depth; i++)
    {
        accum += weight * noise(perlin, temp_p);
        weight *= 0.5f;
        temp_p *= 2.0f;
    }

    return (f32)fabs(accum);
}



#endif
