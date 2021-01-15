
Color emitted(List<Texture>* list, Material* material, const Ray& r, const Hit_Record& record, f32 u, f32 v, const Point3& p)
{
    switch(material->type)
    {
    case MATERIAL_DIFFUSE_LIGHT:
    {
        if (record.front_face)
            return value(list, material->diffuse_light.emit_texture, u, v, p);
        else
            return color(0.0f, 0.0f, 0.0f);
    }
    break;
    default:
    return color(0.0f, 0.0f, 0.0f);
    }
}

f32 scattering_pdf(Material* material, const Ray& r, const Hit_Record& record, Ray& scattered)
{
    switch(material->type)
    {
    case Material_Type::MATERIAL_LAMBERTIAN:
    {
        f32 cosine = dot(record.normal, unit_vector(scattered.direction));
        return cosine < 0.0f ? 0.0f : cosine/pi;
    }
    break;
    }
    return 0.0f;
}

bool scatter(List<Texture>* texture_list, Material* material, const Ray& r, const Hit_Record& record, Scatter_Record& srec)
{
    switch(material->type)
    {
    case Material_Type::MATERIAL_LAMBERTIAN:
    {
        srec.is_specular = false;
        srec.attenuation = value(texture_list, material->lambertian.albedo_handle, record.u, record.v, record.p);
        srec.pdf = alloc_pdf(cosine(record.normal));
        return true;
    }
    case Material_Type::MATERIAL_METAL:
    {
        Vec3 reflected = reflect(unit_vector(r.direction), record.normal);
        srec.specular_ray = ray(record.p, reflected + material->metal.fuzz * random_in_unit_sphere());
        srec.attenuation = material->metal.albedo;
        srec.is_specular = true;
        srec.pdf = nullptr;
        return true;
    }
    case Material_Type::MATERIAL_DIELECTRIC:
    {
        srec.is_specular = true;
        srec.pdf = nullptr;
        srec.attenuation = color(1.0f, 1.0f, 1.0f);
        f32 ir = material->dielectric.ir;
        f32 refraction_ratio = record.front_face ? (1.0f / ir) : ir;

        Vec3 unit_direction = unit_vector(r.direction);
        f32 cos_theta = MIN(dot(-unit_direction, record.normal), 1.0f);
        f32 sin_theta = fsqrt(1.0f - cos_theta * cos_theta);

        bool cannot_refract = refraction_ratio * sin_theta > 1.0f;
        Vec3 direction;

        if(cannot_refract || schlick(cos_theta, refraction_ratio) > random_float())
            direction = reflect(unit_direction, record.normal);
        else
            direction = refract(unit_direction, record.normal, refraction_ratio);
                
        srec.specular_ray = ray(record.p, direction, r.time);
        return true;
    }
    break;
    // case MATERIAL_ISOTROPIC:
    // {
    //     scattered = ray(record.p, random_in_unit_sphere(), r.time);
    //     attenuation = value(texture_list, material->isotropic.albedo_handle, record.u, record.v, record.p);
    //     return true;
    // }
    // break;
    }
    return false;
}
