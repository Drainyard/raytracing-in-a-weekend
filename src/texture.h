#ifndef TEXTURE_H
#define TEXTURE_H

enum Texture_Type
{
    TEXTURE_SOLID_COLOR,
    TEXTURE_CHECKER,
    TEXTURE_NOISE,
    TEXTURE_IMAGE
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
            f32 scale;
        } noise;
        struct
        {
            unsigned char* data;
            i32 width;
            i32 height;
            i32 bytes_per_scanline;
        } image;
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

Texture noise(f32 scale)
{
    Texture t = {};
    t.type = TEXTURE_NOISE;
    t.noise.noise = perlin();
    t.noise.scale = scale;
    return t;
}

const static int BYTES_PER_PIXEL = 3;
Texture image(const char* filename)
{
    i32 components_per_pixel = BYTES_PER_PIXEL;

    Texture t = {};
    t.type = TEXTURE_IMAGE;
    t.image.data = stbi_load(filename, &t.image.width, &t.image.height, &components_per_pixel, components_per_pixel);

    if(!t.image.data)
    {
        printf("ERROR: Could not load texture image file: %s\n", filename);
    }

    t.image.bytes_per_scanline = BYTES_PER_PIXEL * t.image.width;

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
        f32 scale = texture->noise.scale;
        Perlin& noise = texture->noise.noise;
        /* return color(1.0f, 1.0f, 1.0f) * 0.5f * (1.0f + noise(&noise, scale * p)); */
        return color(1.0f, 1.0f, 1.0f) * 0.5f * (1.0f + sin(scale * p.z + 10.0f * turb(&noise, scale * p)));
//        return color(1.0f, 1.0f, 1.0f) * 0.5f * (1.0f + turb(&noise, scale * p));
    }
    break;
    case TEXTURE_IMAGE:
    {
        if(!texture->image.data)
        {
            return color(0.0f, 1.0f, 1.0f);
        }

        u = clamp(u, 0.0f, 1.0f);
        v = 1.0f - clamp(v, 0.0f, 1.0f);

        i32 width = texture->image.width;
        i32 height = texture->image.height;

        i32 i = (i32)(u * width);
        i32 j = (i32)(v * height);

        if(i >= width) i = width - 1;
        if(j >= height) j = height - 1;

        const f32 color_scale = 1.0f / 255.0f;

        unsigned char* data = texture->image.data;

        unsigned char *pixel = data + j * texture->image.bytes_per_scanline + i * BYTES_PER_PIXEL;
        
        return color(color_scale * pixel[0], color_scale * pixel[1], color_scale * pixel[2]);
    }
    break;
    }
    return color(0.0f, 0.0f, 0.0f);
}

#endif
