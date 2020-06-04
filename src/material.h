#ifndef MATERIAL_H
#define MATERIAL_H

enum Material_Type
{
    MATERIAL_LAMBERTIAN,
    MATERIAL_METAL
};

struct Material
{
    Material_Type type;

    union
    {
        struct
        {
            Color albedo;
        } lambertian;
        struct
        {
            Color albedo;
        } metal;
    };
};

#endif
