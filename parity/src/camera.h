/* camera: one pinhole camera with two inverse uses of the same intrinsics.
   camera_ray()      -> path tracer primary rays
   camera_project()  -> rasterizer view-space to pixel coordinates
   Both share eye/basis/fov/aspect, so a point projects to the same pixel
   whichever back end draws it. */
#ifndef PARITY_CAMERA_H
#define PARITY_CAMERA_H

#include "math3.h"

typedef struct {
    V3 eye, forward, right, up;
    float fov_y;       /* vertical fov, radians */
    float aspect;
    int w, h;
} Camera;

Camera camera_make(V3 eye, V3 look_at, V3 up_hint, float fov_y_deg, int w, int h);

/* sx, sy in [0,1) pixel-sample position (jitter lives in the caller) */
Ray camera_ray(const Camera* c, float sx, float sy);

/* project a world point; returns 0 if behind the camera.
   px,py in pixel coordinates (top-left origin), view_z > 0 in front. */
int camera_project(const Camera* c, V3 p, float* px, float* py, float* view_z);

#endif
