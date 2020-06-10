#ifndef TEXTURE_H
#define TEXTURE_H

enum Texture_Type
{
    TEXTURE_SOLID_COLOR,
    TEXTURE_CHECKER,
    TEXTURE_NOISE
};

struct Texture
{
    Texture_Type type;
    union
    {
        struct
        {
            Color color_value;            
        } solid_color;
        struct
        {
            size_t odd_handle;
            size_t even_handle;
        } checker;
        struct
        {
            Perlin noise;
        } noise;
    };
};

Texture solid_color(Color c)
{
    Texture t = {};
    t.type = TEXTURE_SOLID_COLOR;
    t.solid_color.color_value = c;
    return t;
}

Texture solid_color(f32 r, f32 g, f32 b)
{
    return solid_color(color(r, g, b));
}

Texture checkered(size_t t0, size_t t1)
{
    Texture t = {};
    t.type = TEXTURE_CHECKER;
    t.checker.odd_handle = t0;
    t.checker.even_handle = t1;
    return t;
}

Texture noise()
{
    Texture t = {};
    t.type = TEXTURE_NOISE;
    t.noise.noise = perlin();
    return t;
}

Color value(List<Texture>* list, size_t handle, f32 u, f32 v, const Vec3& p)
{
    Texture* texture = &list->data[handle];
    switch(texture->type)
    {
    case TEXTURE_SOLID_COLOR:
    {
        return texture->solid_color.color_value;
    }
    break;
    case TEXTURE_CHECKER:
    {
        f32 sines = fsin(10.0f * p.x) * fsin(10.0f * p.y) * fsin(10.0f * p.z);
        if(sines < 0)
        {
            return value(list, texture->checker.odd_handle, u, v, p);
        }
        else
        {
            return value(list, texture->checker.even_handle, u, v, p);
        }
    }
    break;
    case TEXTURE_NOISE:
    {
        return color(1.0f, 1.0f, 1.0f) * noise(&texture->noise.noise, p);
    }
    break;
    }
    return color(0.0f, 0.0f, 0.0f);
}

#endif
