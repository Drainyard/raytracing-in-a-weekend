
bool scatter(Material* material, const Ray& r, const Hit_Record& record, Color& attenuation, Ray& scattered)
{
    switch(material->type)
    {
    case Material_Type::MATERIAL_LAMBERTIAN:
    {
        Vec3 scatter_direction = record.normal + random_unit_vector();
        scattered.origin = record.p;
        scattered.direction = scatter_direction;
        attenuation = material->lambertian.albedo;
        return true;
    }
    break;
    case Material_Type::MATERIAL_METAL:
    {
        Vec3 reflected = reflect(unit_vector(r.direction), record.normal);
        Ray scatter_ray = {};
        scatter_ray.origin = record.p;
        scatter_ray.direction = reflected;
        scattered = scatter_ray;
        attenuation = material->metal.albedo;
        return (dot(scattered.direction, record.normal) > 0);
    }
    break;
    }
    return false;
}
