#include "camera.h"

Camera camera_make(V3 eye, V3 look_at, V3 up_hint, float fov_y_deg, int w, int h) {
    Camera c;
    c.eye = eye;
    c.forward = vnorm(vsub(look_at, eye));
    c.right = vnorm(vcross(c.forward, up_hint));
    c.up = vcross(c.right, c.forward);
    c.fov_y = fov_y_deg * PI / 180.0f;
    c.aspect = (float)w / (float)h;
    c.w = w; c.h = h;
    return c;
}

Ray camera_ray(const Camera* c, float sx, float sy) {
    float scale_y = tanf(c->fov_y * 0.5f);
    float scale_x = scale_y * c->aspect;
    /* ndc in [-1,1], y flipped so sy=0 is the top row */
    float ndc_x = 2.0f * sx - 1.0f;
    float ndc_y = 1.0f - 2.0f * sy;
    V3 d = vadd(c->forward,
           vadd(vscl(c->right, ndc_x * scale_x),
                vscl(c->up, ndc_y * scale_y)));
    Ray r = { c->eye, vnorm(d) };
    return r;
}

int camera_project(const Camera* c, V3 p, float* px, float* py, float* view_z) {
    V3 rel = vsub(p, c->eye);
    float vz = vdot(rel, c->forward);
    if (vz <= 1e-4f) return 0;
    float rx = vdot(rel, c->right);
    float uy = vdot(rel, c->up);
    float scale_y = tanf(c->fov_y * 0.5f);
    float scale_x = scale_y * c->aspect;
    float ndc_x = (rx / vz) / scale_x;
    float ndc_y = (uy / vz) / scale_y;
    *px = (ndc_x + 1.0f) * 0.5f * c->w;
    *py = (1.0f - ndc_y) * 0.5f * c->h;
    *view_z = vz;
    return 1;
}
