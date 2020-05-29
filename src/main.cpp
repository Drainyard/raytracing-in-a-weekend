#include "types.h"
#include "stdio.h"
#include <stdlib.h>

#include <math.h>
#include "math_util.h"
#include "vec3.h"
#include "color.h"
#include "ray.h"
#include "hittable.h"
#include "camera.h"

Color ray_color(const Ray& r, const Hittable_List* world)
{
    Hit_Record record = {};
    if(hit(world, r, 0.0f, infinity, record))
    {
        return 0.5f * (record.normal + color(1.0f, 1.0f, 1.0f));
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

    FILE* image_file = fopen("output.ppm", "w");

    if(image_file)
    {
        Camera camera = create_camera();
        
        fprintf(image_file, "P3\n%d %d \n255\n", image_width, image_height);

        Hittable_List list = {};
        Hittable h1 = {};
        h1.type = HITTABLE_SPHERE;
        h1.sphere.center = {0.0f, 0.0f, -1.0f};
        h1.sphere.radius = 0.5f;
        add(&list, h1);
        Hittable h2 = {};
        h2.type = HITTABLE_SPHERE;
        h2.sphere.center = {0.0f, -100.5f, -1.0f};
        h2.sphere.radius = 100.0f;
        add(&list, h2);

        for(i32 j = image_height - 1; j >= 0; --j)
        {
            fprintf(stderr, "\rScanlines remaining: %d", j);
            Color pixel_color = {};
            pixel_color.x = 0.0f;
            pixel_color.y = 0.0f;
            pixel_color.z = 0.0f;
            for(i32 i = 0; i < image_width; ++i)
            {
                for(i32 s = 0; s < samples_per_pixel; s++)
                {
                    f32 u = (i + random_float()) / (image_width - 1);
                    f32 v = (j + random_float()) / (image_height - 1);
                    Ray r = get_ray(&camera, u, v);
                    pixel_color += ray_color(r, &list);
                }
                write_color(image_file, pixel_color, samples_per_pixel);
            }
        }
        
        fclose(image_file);
    }
}
