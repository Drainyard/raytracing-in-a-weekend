#include "types.h"
#include "stdio.h"
#include <stdlib.h>
#include <time.h>
#include <assert.h>
#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

#include "list.h"

#include <math.h>
#include "math_util.h"
#include "vec3.h"
#include "onb.h"
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
    f32 pdf;
    Color albedo;

    if(!scatter(texture_list, &material_list->data[record.material_handle], r, record, albedo, scattered, pdf))
    {
        return emitted_color;
    }

    Point3 on_light = point3(random_float(213, 343), 554, random_float(227, 332));
    Vec3 to_light = on_light - record.p;
    f32 distance_squared = length_squared(to_light);
    to_light = unit_vector(to_light);

    if(dot(to_light, record.normal) < 0.0f)
    {
        return emitted_color;
    }

    f32 light_area = (343 - 213) * (332 - 227);
    f32 light_cosine = fabs(to_light.y);
    if(light_cosine < 0.000001f)
    {
        return emitted_color;
    }

    pdf = distance_squared / (light_cosine * light_area);
    scattered = ray(record.p, to_light, r.time);

    return emitted_color +
        albedo * scattering_pdf(&material_list->data[record.material_handle], r, record, scattered)
        * ray_color(texture_list, material_list, scattered, background, world, depth - 1) / pdf;
}

List<Hittable> cornell_box(List<Material>& material_list, List<Texture>& texture_list)
{
    List<Hittable> list = {};

    size_t red = add(&material_list, lambertian(add(&texture_list, solid_color(0.65f, 0.05f, 0.05f))));
    size_t white = add(&material_list, lambertian(add(&texture_list, solid_color(0.73f, 0.73f, 0.73f))));
    size_t green = add(&material_list, lambertian(add(&texture_list, solid_color(0.12f, 0.45f, 0.15f))));
    size_t light = add(&material_list, diffuse_light(add(&texture_list, solid_color(7.0f, 7.0f, 7.0f))));

    add(&list, yz_rect(0, 555, 0, 555, 555, green));
    add(&list, yz_rect(0, 555, 0, 555, 0, red));
    add(&list, xz_rect(213, 343, 227, 332, 554, light));
    add(&list, xz_rect(0, 555, 0, 555, 0, white));
    add(&list, xz_rect(0, 555, 0, 555, 555, white));
    add(&list, xy_rect(0, 555, 0, 555, 555, white));
    
    Hittable* main_box_1 = alloc_hittable(box(point3(0, 0, 0), point3(165, 330, 165), white));

    Hittable* rotate_box = alloc_hittable(rotate_y(&list, main_box_1, 15));
    add(&list, translate(rotate_box, vec3(265, 0, 295)));

    Hittable *main_box_2 = alloc_hittable(box(point3(0, 0, 0), point3(165, 165, 165), white));

    Hittable* rotate_box_2 = alloc_hittable(rotate_y(&list, main_box_2, -18));
    add(&list, translate(rotate_box_2, vec3(130, 0, 65)));
    
    
    return list;
}

List<Hittable> cornell_smoke(List<Material>& material_list, List<Texture>& texture_list)
{
    List<Hittable> list = {};

    size_t red = add(&material_list, lambertian(add(&texture_list, solid_color(0.65f, 0.05f, 0.05f))));
    size_t white = add(&material_list, lambertian(add(&texture_list, solid_color(0.73f, 0.73f, 0.73f))));
    size_t green = add(&material_list, lambertian(add(&texture_list, solid_color(0.12f, 0.45f, 0.15f))));
    size_t light = add(&material_list, diffuse_light(add(&texture_list, solid_color(7.0f, 7.0f, 7.0f))));

    add(&list, yz_rect(0, 555, 0, 555, 555, green));
    add(&list, yz_rect(0, 555, 0, 555, 0, red));
    add(&list, xz_rect(113, 443, 127, 432, 554, light));
    add(&list, xz_rect(0, 555, 0, 555, 0, white));
    add(&list, xz_rect(0, 555, 0, 555, 555, white));
    add(&list, xy_rect(0, 555, 0, 555, 555, white));
    
    Hittable* main_box_1 = alloc_hittable(box(point3(0, 0, 0), point3(165, 330, 165), white));

    Hittable* rotate_box = alloc_hittable(rotate_y(&list, main_box_1, 15));
    Hittable* t_1 = alloc_hittable(translate(rotate_box, vec3(265, 0, 295)));

    Hittable *main_box_2 = alloc_hittable(box(point3(0, 0, 0), point3(165, 165, 165), white));

    Hittable* rotate_box_2 = alloc_hittable(rotate_y(&list, main_box_2, -18));
    Hittable* t_2 = alloc_hittable(translate(rotate_box_2, vec3(130, 0, 65)));

    add(&list, constant_medium(&material_list, t_1, 0.01f, add(&texture_list, solid_color(color(0, 0, 0)))));
    add(&list, constant_medium(&material_list, t_2, 0.01f, add(&texture_list, solid_color(color(1, 1, 1)))));
    
    return list;
}

List<Hittable> final_scene(List<Material>& material_list, List<Texture>& texture_list)
{
    List<Hittable> boxes1 = {};

    size_t ground = add(&material_list, lambertian(add(&texture_list, solid_color(0.48, 0.83, 0.53))));

    const i32 boxes_per_side = 20;
    for(i32 i = 0; i < boxes_per_side; i++)
    {
        for(i32 j = 0; j < boxes_per_side; j++)
        {
            f32 w = 100.0f;
            f32 x0 = -1000.0f + i * w;
            f32 z0 = -1000.0f + j * w;
            f32 y0 = 0.0f;
            f32 x1 = x0 + w;
            f32 y1 = random_float(1, 101);
            f32 z1 = z0 + w;

            add(&boxes1, box(point3(x0, y0, z0), point3(x1, y1, z1), ground));
        }
    }   

    List<Hittable> objects = {};

    add(&objects, bvh_node(&boxes1, 0, 1));

    size_t light = add(&material_list, diffuse_light(add(&texture_list, solid_color(7, 7, 7))));
    add(&objects, xz_rect(123, 423, 147, 412, 554, light));

    Point3 center1 = point3(400, 400, 200);
    Point3 center2 = center1 + vec3(30, 0, 0);
    size_t moving_sphere_material = add(&material_list, lambertian(add(&texture_list, solid_color(0.7, 0.3, 0.1))));

    add(&objects, moving_sphere(center1, center2, 0, 1, 50, moving_sphere_material));

    add(&objects, sphere(point3(260, 150, 45), 50, add(&material_list, dialectric(1.5f))));
    add(&objects, sphere(point3(0, 150, 145), 50, add(&material_list, metal(color(0.8f, 0.8f, 0.9f), 1.0f))));

    size_t boundary = add(&objects, sphere(point3(360, 150, 145), 70, add(&material_list, dialectric(1.5f))));
    add(&objects, constant_medium(&material_list, &objects.data[boundary], 0.2f, add(&texture_list, solid_color(0.2f, 0.4f, 0.9f))));

    Hittable* sphere_boundary = alloc_hittable(sphere(point3(0, 0, 0), 5000, add(&material_list, dialectric(1.5f))));
    add(&objects, constant_medium(&material_list, sphere_boundary, 0.0001f, add(&texture_list, solid_color(1, 1, 1))));

    size_t emat = add(&material_list, lambertian(add(&texture_list, image("earthmap.jpg"))));
    add(&objects, sphere(point3(400, 200, 400), 100, emat));

    size_t pertext = add(&texture_list, noise(0.1f));
    add(&objects, sphere(point3(220, 280, 300), 80, add(&material_list, lambertian(pertext))));
    
    List<Hittable> boxes2 = {};

    size_t white = add(&material_list, lambertian(add(&texture_list, solid_color(.73f, .73f, .73f))));
    i32 ns = 1000;
    for(i32 j = 0; j < ns; j++)
    {
        add(&boxes2, sphere(random_vec3(0, 165), 10, white));
    }

    Hittable* bvh = alloc_hittable(bvh_node(&boxes2, 0.0f, 1.0f));
    Hittable* rotated = alloc_hittable(rotate_y(&objects, bvh, 15.0f));
    add(&objects, translate(rotated, vec3(-100, 290, 395)));

    return objects;
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

i32 main()
{
    time_t t;
    srand((unsigned) time(&t));

    f32 aspect_ratio = 1.0f / 1.0f;
    i32 image_width = 600;
    i32 image_height = i32(image_width / aspect_ratio);
    const i32 samples_per_pixel = 100;
    const i32 max_depth = 50;

    char filename[256];
    sprintf(filename, "%d_samples_cornell.ppm", samples_per_pixel);

    FILE* image_file = fopen(filename, "w");

    if(image_file)
    {
        Vec3 look_from = point3(278, 278, -800);
        Vec3 look_at = point3(278, 278, 0.0f);
        Vec3 v_up = vec3(0.0f, 1.0f, 0.0f);
        f32 dist_to_focus = 10.0f;
        f32 aperture = 0.0f;
        f32 v_fov = 40.0f;       

        List<Material> material_list = {};
        List<Texture> texture_list = {};
        List<Hittable> list = {};

        Color background = color(0.0f, 0.0f, 0.0f);

        switch(6)
        {
        case 1:
        {
            aspect_ratio = 16.0f / 9.0f;
            image_width = 600;
            image_height = i32(image_width / aspect_ratio);
            
            list = random_scene(material_list, texture_list);
            look_from = point3(13,2,3);
            look_at = point3(0,0,0);
            v_fov = 20.0;
            background = color(0.70f, 0.80f, 1.00f);
        }
        break;

        case 2:
        {
            aspect_ratio = 16.0f / 9.0f;
            image_width = 600;
            image_height = i32(image_width / aspect_ratio);
            
            list = two_spheres(material_list, texture_list);
            look_from = point3(13,2,3);
            look_at = point3(0,0,0);
            v_fov = 20.0;
            background = color(0.70f, 0.80f, 1.00f);
        }
        break;

        case 3:
        {
            aspect_ratio = 16.0f / 9.0f;
            image_width = 600;
            image_height = i32(image_width / aspect_ratio);
            
            list = two_perlin_spheres(material_list, texture_list);
            look_from = point3(13,2,3);
            look_at = point3(0,0,0);
            v_fov = 20.0;
            background = color(0.70f, 0.80f, 1.00f);
        }
        break;

        case 4:
        {
            aspect_ratio = 16.0f / 9.0f;
            image_width = 600;
            image_height = i32(image_width / aspect_ratio);

            look_from = point3(0, 0, 12);
            look_at = point3(0, 0, 0);
            v_fov = 20.0f;
            list = earth(material_list, texture_list);    
        }
        break;
        case 5:
        {
            aspect_ratio = 16.0f / 9.0f;
            image_width = 600;
            image_height = i32(image_width / aspect_ratio);
            
            look_from = point3(26, 3, 6);
            look_at = point3(0, 2, 0);
            v_fov = 20.0f;
            list = simple_light(material_list, texture_list);    
        }
        break;
        case 6:
        {
            aspect_ratio = 1.0f;
            image_width = 600;
            image_height = i32(image_width / aspect_ratio);
            
            look_from = point3(278, 278, -800);
            look_at = point3(278, 278, 0);
            v_fov = 40.0f;
            list = cornell_box(material_list, texture_list);
        }
        break;
        case 7:
        {
            aspect_ratio = 1.0f;
            image_width = 400;
            image_height = i32(image_width / aspect_ratio);
            
            look_from = point3(278, 278, -800);
            look_at = point3(278, 278, 0);
            v_fov = 40.0f;
            list = cornell_smoke(material_list, texture_list);
        }
        break;
        case 8:
        {
            aspect_ratio = 1.0f;
            image_width = 600;
            image_height = i32(image_width / aspect_ratio);
            
            look_from = point3(478, 278, -600);
            look_at = point3(278, 278, 0);
            v_fov = 40.0f;
            list = final_scene(material_list, texture_list);
        }
        break;
        }

        fprintf(image_file, "P3\n%d %d \n255\n", image_width, image_height);
        Camera camera = create_camera(look_from, look_at, v_up,
                                      v_fov, aspect_ratio, aperture, dist_to_focus, 0.0f, 1.0f);

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
        fprintf(stderr, "\n");
        
        fclose(image_file);
    }
}
