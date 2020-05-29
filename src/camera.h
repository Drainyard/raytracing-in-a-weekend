#ifndef CAMERA_H
#define CAMERA_H

struct Camera
{
    Point3 origin;
    Point3 lower_left_corner;
    Vec3 horizontal;
    Vec3 vertical;
};

Camera create_camera()
{
    Camera camera = {};

    f32 aspect_ratio = 16.0f / 9.0f;
    f32 viewport_height = 2.0f;
    f32 viewport_width = aspect_ratio * viewport_height;
    f32 focal_length = 1.0f;

    camera.origin = {0.0f, 0.0f, 0.0f};
    camera.horizontal = vec3(viewport_width, 0.0f, 0.0f);
    camera.vertical = vec3(0.0f, viewport_height, 0.0f);
    camera.lower_left_corner = camera.origin - camera.horizontal/2 - camera.vertical/2 - vec3(0.0f, 0.0f, focal_length);
    return camera;
}

Ray get_ray(Camera* camera, f32 u, f32 v)
{
    Ray ray = {};
    ray.origin = camera->origin;
    ray.direction = camera->lower_left_corner + u * camera->horizontal + v * camera->vertical - camera->origin;
    return ray;
}

#endif
