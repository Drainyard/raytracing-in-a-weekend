
Color emitted(List<Texture>* list, Material* material, f32 u, f32 v, const Point3& p)
{
    switch(material->type)
    {
    case MATERIAL_DIFFUSE_LIGHT:
    {
        return value(list, material->diffuse_light.emit_texture, u, v, p);
    }
    break;
    default:
    return color(0.0f, 0.0f, 0.0f);
    }
}

bool scatter(List<Texture>* texture_list, Material* material, const Ray& r, const Hit_Record& record, Color& attenuation, Ray& scattered)
{
    switch(material->type)
    {
    case Material_Type::MATERIAL_LAMBERTIAN:
    {
        Vec3 scatter_direction = record.normal + random_unit_vector();
        scattered = ray(record.p, scatter_direction, r.time);
        attenuation = value(texture_list, material->lambertian.albedo_handle, record.u, record.v, record.p);
        return true;
    }
    break;
    case Material_Type::MATERIAL_METAL:
    {
        Vec3 reflected = reflect(unit_vector(r.direction), record.normal);
        scattered = ray(record.p, reflected + material->metal.fuzz * random_in_unit_sphere());
        attenuation = material->metal.albedo;
        return (dot(scattered.direction, record.normal) > 0);
    }
    break;
    case Material_Type::MATERIAL_DIALECTRIC:
    {
        attenuation = color(1.0f, 1.0f, 1.0f);
        f32 etai_over_etat;
        if(record.front_face)
        {
            etai_over_etat = 1.0f / material->dialectric.ref_idx;
        }
        else
        {
            etai_over_etat = material->dialectric.ref_idx;
        }

        Vec3 unit_direction = unit_vector(r.direction);
        f32 cos_theta = MIN(dot(-unit_direction, record.normal), 1.0f);
        f32 sin_theta = fsqrt(1.0f - cos_theta * cos_theta);
        if(etai_over_etat * sin_theta > 1.0f)
        {
            Vec3 reflected = reflect(unit_direction, record.normal);
            scattered = ray(record.p, reflected);
            return true;
        }

        f32 reflect_prob = schlick(cos_theta, etai_over_etat);
        if(random_float() < reflect_prob)
        {
            Vec3 reflected = reflect(unit_direction, record.normal);
            scattered = ray(record.p, reflected);
            return true;
        }
        
        Vec3 refracted = refract(unit_direction, record.normal, etai_over_etat);
        scattered = ray(record.p, refracted);
        return true;
    }
    break;
    case MATERIAL_ISOTROPIC:
    {
        scattered = ray(record.p, random_in_unit_sphere(), r.time);
        attenuation = value(texture_list, material->isotropic.albedo_handle, record.u, record.v, record.p);
        return true;
    }
    break;
    }
    return false;
}
