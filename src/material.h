#ifndef MATERIAL_H
#define MATERIAL_H

enum Material_Type
{
    MATERIAL_LAMBERTIAN,
    MATERIAL_METAL,
    MATERIAL_DIALECTRIC,
    MATERIAL_DIFFUSE_LIGHT,
    MATERIAL_ISOTROPIC
};

struct Material
{
    Material_Type type;

    union
    {
        struct
        {
            size_t albedo_handle;
        } lambertian;
        struct
        {
            Color albedo;
            f32 fuzz;
        } metal;
        struct
        {
            f32 ref_idx;
        } dialectric;
        struct
        {
            size_t emit_texture;
        } diffuse_light;
        struct
        {
            size_t albedo_handle;
        } isotropic;
    };
};

f32 schlick(f32 cosine, f32 ref_idx)
{
    f32 r0 = (1.0f - ref_idx) / (1.0f + ref_idx);
    r0 = r0 * r0;
    return r0 + (1.0f - r0) * fpow((1.0f - cosine), 5);
}

Material lambertian(size_t albedo)
{
    Material material = {};
    material.type = MATERIAL_LAMBERTIAN;
    material.lambertian.albedo_handle = albedo;
    return material;
}

Material metal(Color albedo, f32 fuzz)
{
    Material material = {};
    material.type = MATERIAL_METAL;
    material.metal.albedo = albedo;
    material.metal.fuzz = fuzz;
    return material;
}

Material dialectric(f32 ref_idx)
{
    Material material = {};
    material.type = MATERIAL_DIALECTRIC;
    material.dialectric.ref_idx = ref_idx;
    return material;
}

Material diffuse_light(size_t emit)
{
    Material material = {};
    material.type = MATERIAL_DIFFUSE_LIGHT;
    material.diffuse_light.emit_texture = emit;
    return material;
}

Material isotropic(size_t texture)
{
    Material material = {};
    material.type = MATERIAL_ISOTROPIC;
    material.isotropic.albedo_handle = texture;
    return material;
}

#endif
