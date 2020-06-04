#include "types.h"
#include "stdio.h"
#include <stdlib.h>

#include <math.h>
#include "math_util.h"
#include "vec3.h"
#include "color.h"
#include "ray.h"
#include "material.h"
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

int main()
{
    const f32 aspect_ratio = 16.0f / 9.0f;
    const int image_width = 384 * 2;
    const int image_height = i32(image_width / aspect_ratio);
    const int samples_per_pixel = 100;
    const int max_depth = 50;

    FILE* image_file = fopen("output.ppm", "w");

    if(image_file)
    {
        Camera camera = create_camera();
        
        fprintf(image_file, "P3\n%d %d \n255\n", image_width, image_height);

        Hittable_List list = {};

        Material mat_1 = {};
        mat_1.lambertian.albedo = color(0.7f, 0.3f, 0.3f);
        mat_1.type = MATERIAL_LAMBERTIAN;
        Material mat_2 = {};
        mat_2.lambertian.albedo = color(0.8f, 0.8f, 0.0f);
        mat_2.type = MATERIAL_LAMBERTIAN;

        Material mat_3 = {};
        mat_3.metal.albedo = color(0.6f, 0.8f, 0.2f);
        mat_3.type = MATERIAL_METAL;
        Material mat_4 = {};
        mat_4.lambertian.albedo = color(0.8f, 0.8f, 0.8f);
        mat_4.type = MATERIAL_METAL;        
        

        Hittable h1 = {};
        h1.type = HITTABLE_SPHERE;
        h1.sphere.center = {0.0f, 0.0f, -1.0f};
        h1.sphere.radius = 0.5f;
        h1.material = &mat_1;
        add(&list, h1);

        Hittable h2 = {};
        h2.type = HITTABLE_SPHERE;
        h2.sphere.center = {0.0f, -100.5f, -1.0f};
        h2.sphere.radius = 100.0f;
        h2.material = &mat_2;
        add(&list, h2);

        Hittable h3 = {};
        h3.type = HITTABLE_SPHERE;
        h3.sphere.center = {1.0f, 0.0f, -1.0f};
        h3.sphere.radius = 0.5f;
        h3.material = &mat_3;
        add(&list, h3);

        Hittable h4 = {};
        h4.type = HITTABLE_SPHERE;
        h4.sphere.center = {-1.0f, 0.0f, -1.0f};
        h4.sphere.radius = 0.5f;
        h4.material = &mat_4;
        add(&list, h4);


        

        for(i32 j = image_height - 1; j >= 0; --j)
        {
            fprintf(stderr, "\rScanlines remaining: %d", j);
            for(i32 i = 0; i < image_width; ++i)
            {
                Color pixel_color = {};
                pixel_color.x = 0.0f;
                pixel_color.y = 0.0f;
                pixel_color.z = 0.0f;
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
