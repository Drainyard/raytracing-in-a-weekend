#include "types.h"
#include "stdio.h"
#include <stdlib.h>
#include <assert.h>
#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

#include "list.h"

#include <math.h>
#include "math_util.h"
#include "vec3.h"
#include "color.h"
#include "ray.h"
#include "perlin.h"
#include "texture.h"
#include "material.h"
#include "aabb.h"
#include "hittable.h"
#include "camera.h"
#include "material.cpp"

Color ray_color(List<Texture>* texture_list, List<Material>* material_list, const Ray& r, const Color background, const List<Hittable>* world, i32 depth)
{
    Hit_Record record = {};
    if(depth <= 0)
    {
        return color(0.0f, 0.0f, 0.0f);
                                           
    }

    if(!hit(world, r, 0.001f, infinity, record))
    {
        return background;
    }
    

    Ray scattered = {};
    Color attenuation = {};
    Color emitted_color = emitted(texture_list, &material_list->data[record.material_handle], record.u, record.v, record.p);

    if(!scatter(texture_list, &material_list->data[record.material_handle], r, record, attenuation, scattered))
    {
        return emitted_color;
    }
        
    return emitted_color + attenuation * ray_color(texture_list, material_list, scattered, background, world, depth - 1);
}

List<Hittable> simple_light(List<Material>& material_list, List<Texture>& texture_list)
{
    List<Hittable> list = {};

    size_t per_text = add(&texture_list, noise(10.0f));
    add(&list, sphere(point3(0.0f, -1000.0f, 0.0f), 1000.0f, add(&material_list, lambertian(per_text))));
    add(&list, sphere(point3(0.0f, 2.0f, 0.0f), 2.0f, add(&material_list, lambertian(per_text))));

    size_t diff_light = add(&material_list, diffuse_light(add(&texture_list, solid_color(4.0f, 4.0f, 4.0f))));
    add(&list, sphere(point3(0.0f, 7.0f, 0.0f), 2.0f, diff_light));
    add(&list, xy_rect(3, 5, 1, 3, -2, diff_light));

    return list;
}

List<Hittable> two_spheres(List<Material>& material_list, List<Texture>& texture_list)
{
    List<Hittable> list = {};

    size_t s0 = add(&texture_list, solid_color(color(0.2f, 0.3f, 0.1f)));
    size_t s1 = add(&texture_list, solid_color(color(0.9f, 0.9f, 0.9f)));

    size_t checker = add(&texture_list, checkered(s0, s1));

    add(&list, sphere(point3(0.0f, -10.0f, 0.0f), 10.0f, add(&material_list, lambertian(checker))));
    add(&list, sphere(point3(0.0f, 10.0f, 0.0f), 10.0f, add(&material_list, lambertian(checker))));

    return list;
}

List<Hittable> two_perlin_spheres(List<Material>& material_list, List<Texture>& texture_list)
{
    List<Hittable> list = {};

    size_t per_text = add(&texture_list, noise(4.0f));
    add(&list, sphere(point3(0.0f, -1000.0f, 0.0f), 1000.0f, add(&material_list, lambertian(per_text))));
    add(&list, sphere(point3(0.0f, 2.0f, 0.0f), 2.0f, add(&material_list, lambertian(per_text))));

    return list;
}

List<Hittable> earth(List<Material>& material_list, List<Texture>& texture_list)
{
    List<Hittable> list = {};

    size_t earth_texture = add(&texture_list, image("earthmap.jpg"));
    size_t earth_surface = add(&material_list, lambertian(earth_texture));

    add(&list, sphere(point3(0.0f, 0.0f, 0.0f), 2.0f, earth_surface));

    return list;
}

List<Hittable> random_scene(List<Material>& material_list, List<Texture>& texture_list)
{
    List<Hittable> list = {};

    size_t t0 = add(&texture_list, solid_color(color(0.2f, 0.3f, 0.1f)));
    size_t t1 = add(&texture_list, solid_color(color(0.9f, 0.9f, 0.9f)));
    size_t checker = add(&texture_list, checkered(t0, t1));
    
    Material ground_material = lambertian(checker);
    Hittable ground = sphere(point3(0, -1000.0f, 0.0f), 1000.0f, add(&material_list, ground_material));
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
                    Material mat = lambertian(add(&texture_list, solid_color(albedo)));

                    Point3 center2 = center + vec3(0.0f, random_float(0.0f, 0.5f), 0.0f);
                    
                    Hittable s = moving_sphere(center, center2, 0.0f, 1.0f, 0.2f, add(&material_list, mat));
                    add(&list, s);
                }
                else if(choose_mat < 0.95f)
                {
                    Color albedo = random_vec3(0.5f, 1.0f);
                    f32 fuzz = random_float(0.0f, 0.5f);
                    Material mat = metal(albedo, fuzz);
                    Hittable s = sphere(center, 0.2f, add(&material_list, mat));
                    add(&list, s);
                }
                else
                {
                    Material mat = dialectric(1.5f);
                    Hittable s = sphere(center, 0.2f, add(&material_list, mat));
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
    const int image_width = 600;
    const int image_height = i32(image_width / aspect_ratio);
    const int samples_per_pixel = 1000;
    const int max_depth = 500;

    FILE* image_file = fopen("output.ppm", "w");

    if(image_file)
    {
        Point3 look_from = point3(25.0f, 8.0f, 3.0f);
        Point3 look_at = point3(0.0f, 0.0f, 0.0f);
        Vec3 v_up = vec3(0.0f, 1.0f, 0.0f);
        f32 dist_to_focus = 10.0f;
        f32 aperture = 0.0f;
        Camera camera = create_camera(look_from, look_at, v_up,
                                      20.0f, aspect_ratio, aperture, dist_to_focus, 0.0f, 1.0f);
        
        fprintf(image_file, "P3\n%d %d \n255\n", image_width, image_height);

        List<Material> material_list = {};
        List<Texture> texture_list = {};
        //List<Hittable> list = random_scene(material_list, texture_list);
        //List<Hittable> list = earth(material_list, texture_list);
        List<Hittable> list = simple_light(material_list, texture_list);
        //List<Hittable> list = two_perlin_spheres(material_list, texture_list);

        // Hittable h1 = sphere({0.0f, 1.0f, 0.0f}, 1.0f, add(&material_list, dialectric(1.5f)));
        // add(&list, h1);

        // size_t t2 = add(&texture_list, solid_color(color(0.4f, 0.2f, 0.1f)));
        // Hittable h2 = sphere({-4.0f, 1.0f, 0.0f}, 1.0f, add(&material_list, lambertian(t2)));
        // add(&list, h2);

        // Hittable h3 = sphere({4.0f, 1.0f, 0.0f}, 1.0f, add(&material_list, metal(color(0.7f, 0.6, 0.5f), 0.0f)));
        // add(&list, h3);

        Color background = color(0.0f, 0.0f, 0.0f);

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
                    pixel_color += ray_color(&texture_list, &material_list, r, background, &list, max_depth);
                }
                write_color(image_file, pixel_color, samples_per_pixel);
            }
        }
        
        fclose(image_file);
    }
}
