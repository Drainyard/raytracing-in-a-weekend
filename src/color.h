#ifndef COLOR_H
#define COLOR_H

void write_color(FILE* file, Color pixel_color)
{
    fprintf(file, "%d %d %d\n", i32(255.99 * pixel_color.x),
           i32(255.99 * pixel_color.y),
           i32(255.99 * pixel_color.z));
}

#endif
