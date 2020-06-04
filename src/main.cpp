#include "types.h"
#include "stdio.h"
#include <stdlib.h>

#include <math.h>
#include "math_util.h"
#include "vec3.h"
#include "color.h"
#include "ray.h"
#include "material.h"
#include "list.h"
#include "hittable.h"
#include "camera.h"
#include "material.cpp"

Color ray_color(const Ray& r, const Hittable_List* world, i32 depth)
{
    Hit_Record record = {};
    if(depth <= 0)
    {
        return color(0.0f, 0.0f, 0.0f);
                                           
    }
    if(hit(world, r, 0.001f, infinity, record))
    {
        Ray scattered = {};
        Color attenuation = {};
        if(scatter(record.material, r, record, attenuation, scattered))
        {
            return attenuation * ray_color(scattered, world, depth - 1);
        }
        return color(0.0f, 0.0f, 0.0f);
    }
    Vec3 unit_direction = unit_vector(r.direction);
    f32 t = 0.5f * (unit_direction.y + 1.0f);
    return (1.0f - t) * color(1.0f, 1.0f, 1.0f) + t * color(0.5f, 0.7f, 1.0f);
}

Hittable_List random_scene(List<Material>& material_list)
{
    Hittable_List list = {};

    Material ground_material = lambertian(color(0.5f, 0.5f, 0.5f));
    Material& ground_mat = add(&material_list, ground_material);
    Hittable ground = sphere(point3(0, -1000.0f, 0.0f), 1000.0f, &ground_mat);
    add(&list, ground);

    for(i32 a = -11; a < 11; a++)
    {
        for(i32 b = -11; b < 11; b++)
        {
            f32 choose_mat = random_float();
            Point3 center = point3(a + 0.9f * random_float(), 0.2f, b + 0.9f * random_float());

            if(length(center - vec3(4.0f, 0.2f, 0.0f)) > 0.9f)
            {
                if(choose_mat < 0.8f)
                {
                    Color albedo = random_vec3() * random_vec3();
                    Material mat = lambertian(albedo);
                    Material& m = add(&material_list, mat);
                    Hittable s = sphere(center, 0.2f, &m);
                    add(&list, s);
                }
                else if(choose_mat < 0.95f)
                {
                    Color albedo = random_vec3(0.5f, 1.0f);
                    f32 fuzz = random_float(0.0f, 0.5f);
                    Material mat = metal(albedo, fuzz);
                    Material& m = add(&material_list, mat);
                    Hittable s = sphere(center, 0.2f, &m);
                    add(&list, s);
                }
                else
                {
                    Material mat = dialectric(1.5f);
                    Material& m = add(&material_list, mat);
                    Hittable s = sphere(center, 0.2f, &m);
                    add(&list, s);
                }
            }
        }
    }
    return list;
}

int main()
{
    const f32 aspect_ratio = 16.0f / 9.0f;
    const int image_width = 1200.0f;
    const int image_height = i32(image_width / aspect_ratio);
    const int samples_per_pixel = 100;
    const int max_depth = 50;

    FILE* image_file = fopen("output.ppm", "w");

    if(image_file)
    {
        Point3 look_from = point3(13.0f, 2.0f, 3.0f);
        Point3 look_at = point3(0.0f, 0.0f, 0.0f);
        Vec3 v_up = vec3(0.0f, 1.0f, 0.0f);
        f32 dist_to_focus = 10.0f;
        f32 aperture = 0.1f;
        Camera camera = create_camera(look_from, look_at, v_up,
                                      20.0f, aspect_ratio, aperture, dist_to_focus);
        
        fprintf(image_file, "P3\n%d %d \n255\n", image_width, image_height);

        List<Material> material_list = {};
        Hittable_List list = random_scene(material_list);

        Material m1 = dialectric(1.5f);
        Hittable h1 = sphere({0.0f, 1.0f, 0.0f}, 1.0f, &m1);
        add(&list, h1);

        Material m2 = lambertian(color(0.4f, 0.2f, 0.1f));
        Hittable h2 = sphere({-4.0f, 1.0f, 0.0f}, 1.0f, &m2);
        add(&list, h2);

        Material m3 = metal(color(0.7f, 0.6, 0.5f), 0.0f);
        Hittable h3 = sphere({4.0f, 1.0f, 0.0f}, 1.0f, &m3);
        add(&list, h3);

        for(i32 j = image_height - 1; j >= 0; --j)
        {
            fprintf(stderr, "\rScanlines remaining: %d", j);
            for(i32 i = 0; i < image_width; ++i)
            {
                Color pixel_color = color(0.0f, 0.0f, 0.0f);
                for(i32 s = 0; s < samples_per_pixel; s++)
                {
                    f32 u = (i + random_float()) / (image_width - 1);
                    f32 v = (j + random_float()) / (image_height - 1);
                    Ray r = get_ray(&camera, u, v);
                    pixel_color += ray_color(r, &list, max_depth);
                }
                write_color(image_file, pixel_color, samples_per_pixel);
            }
        }
        
        fclose(image_file);
    }
}
