#ifndef PERLIN_H
#define PERLIN_H

static const int POINT_COUNT = 256;
struct Perlin
{
    f32* rand_float;
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
    perlin.rand_float = (f32*)malloc(sizeof(f32) * POINT_COUNT);
    for(i32 i = 0; i < POINT_COUNT; i++)
    {
        perlin.rand_float[i] = random_float();
    }

    perlin.perm_x = perlin_generate_perm();
    perlin.perm_y = perlin_generate_perm();
    perlin.perm_z = perlin_generate_perm();
    
    return perlin;                  
}

f32 noise(Perlin* perlin, const Point3& p)
{
    f32 u = p.x - floor(p.x);
    f32 v = p.y - floor(p.y);
    f32 w = p.z - floor(p.z);

    i32 i = (i32)(4 * p.x) & 255;
    i32 j = (i32)(4 * p.y) & 255;
    i32 k = (i32)(4 * p.z) & 255;

    return perlin->rand_float[perlin->perm_x[i] ^ perlin->perm_y[j] ^ perlin->perm_z[k]];
}



#endif
