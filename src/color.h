#ifndef COLOR_H
#define COLOR_H

void write_color(FILE* file, Color pixel_color, i32 samples_per_pixel)
{
    f32 r = pixel_color.x;
    f32 g = pixel_color.y;
    f32 b = pixel_color.z;

    f32 scale = 1.0f / samples_per_pixel;
    r *= scale;
    g *= scale;
    b *= scale;
    
    fprintf(file, "%d %d %d\n",
            i32(256 * clamp(r, 0.0f, 0.999f)),
            i32(256 * clamp(g, 0.0f, 0.999f)),
            i32(256 * clamp(b, 0.0f, 0.999f)));
}

#endif
