#include "types.h"
#include "stdio.h"

#include <math.h>
#include "vec3.h"
#include "color.h"
#include "ray.h"

Color ray_color(const Ray& r)
{
    Vec3 unit_direction = unit_vector(r.direction);
    f32 t = 0.5f * (unit_direction.y + 1.0f);
    return (1.0f - t) * color(1.0f, 1.0f, 1.0f) + t * color(0.5f, 0.7f, 1.0f);
}


int main()
{
    const f32 aspect_ratio = 16.0f / 9.0f;
    const int image_width = 384;
    const int image_height = i32(image_width / aspect_ratio);

    FILE* image_file = fopen("output.ppm", "w");

    if(image_file)
    {
        f32 viewport_height = 2.0f;
        f32 viewport_width = aspect_ratio * viewport_height;
        f32 focal_length = 1.0f;

        Point3 origin = {0.0f, 0.0f, 0.0f};
        Vec3 horizontal = vec3(viewport_width, 0.0f, 0.0f);
        Vec3 vertical = vec3(0.0f, viewport_height, 0.0f);
        Vec3 lower_left_corner = origin - horizontal / 2 - vertical / 2 - vec3(0, 0, focal_length);
        
        fprintf(image_file, "P3\n %d %d \n255\n", image_width, image_height);

        for(i32 j = image_height - 1; j >= 0; --j)
        {
            fprintf(stderr, "\rScanlines remaining: %d", j);
            for(i32 i = 0; i < image_width; ++i)
            {
                f32 u = f32(i) / (image_width -1);
                f32 v = f32(j) / (image_height -1);
                Ray r = {};
                r.origin = origin;
                r.direction = lower_left_corner + u * horizontal + v * vertical - origin;
                Color pixel_color = ray_color(r);
                write_color(image_file, pixel_color);
            }
        }
        
        fclose(image_file);
    }
}
