#ifndef CAMERA_H
#define CAMERA_H

struct Camera
{
    Point3 origin;
    Point3 lower_left_corner;
    Vec3 horizontal;
    Vec3 vertical;
    f32 lens_radius;
    Vec3 u, v, w;
    f32 time0;
    f32 time1;
};

Camera create_camera(Point3 look_from, Point3 look_at, Vec3 vup, f32 vfov, f32 aspect_ratio, f32 aperture, f32 focus_dist, f32 t0 = 0.0f, f32 t1 = 0.0f)
{
    Camera camera = {};

    f32 theta = degrees_to_radians(vfov);
    f32 h = ftan(theta * 0.5f);
    
    f32 viewport_height = 2.0f * h;
    f32 viewport_width = aspect_ratio * viewport_height;
    
    camera.w = unit_vector(look_from - look_at);
    camera.u = unit_vector(cross(vup, camera.w));
    camera.v = cross(camera.w, camera.u);

    camera.origin = look_from;
    camera.horizontal = focus_dist * viewport_width * camera.u;
    camera.vertical = focus_dist * viewport_height * camera.v;
    camera.lower_left_corner = camera.origin - camera.horizontal/2 - camera.vertical/2 - focus_dist * camera.w;
    camera.lens_radius = aperture * 0.5f;
    camera.time0 = t0;
    camera.time1 = t1;
    return camera;
}

Ray get_ray(Camera* camera, f32 s, f32 t)
{
    Vec3 rd = camera->lens_radius * random_in_unit_disk();
    Vec3 offset = camera->u * rd.x + camera->v * rd.y;
    return ray(camera->origin + offset,
               camera->lower_left_corner + s * camera->horizontal + t * camera->vertical - camera->origin - offset,
               random_float(camera->time0, camera->time1));
}

#endif
