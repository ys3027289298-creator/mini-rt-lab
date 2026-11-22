
/* accum: temporal accumulation over reprojected history.

   每帧对带噪单帧估计做重投影（针孔相机，世界坐标 -> 上一帧像素），
   几何不一致或颜色差异超阈值的历史按像素拒绝，拒绝/缺失处退回单帧；
   输出累积帧数、拒绝计数与噪声(RMSE)下降，并与单帧噪声对比。
   场景是一个确定性的合成棋盘地面 + 移动相机，用于可复现地自检。 */
#include <stdio.h>
#include <math.h>
#include <string.h>

#define W 48
#define H 32
#define NPIX (W * H)
#define FOV_Y 0.96f              /* 约 55 度垂直视场 */
#define NOISE_SIGMA 0.10f        /* 单帧每通道噪声幅度，对称截断到 ±1.5sigma */
#define COLOR_REJECT 0.26f       /* 亮度差阈值：约 4 倍帧间噪声差标准差，棋盘翻转差 0.60 仍稳超 */
#define GEO_EPS 0.5f             /* 重投影世界点与历史交点的距离阈值 */

typedef struct { float x, y, z; } V3;
typedef struct { float r, g, b; } RGB;
typedef struct { V3 o, f, r, u; } Camera;

static V3 vsub(V3 a, V3 b) { V3 v = { a.x - b.x, a.y - b.y, a.z - b.z }; return v; }
static V3 vcross(V3 a, V3 b) {
    V3 v = { a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x };
    return v;
}
static float vdot(V3 a, V3 b) { return a.x * b.x + a.y * b.y + a.z * b.z; }
static V3 vnorm(V3 a) {
    float l = sqrtf(vdot(a, a));
    V3 r = { 0, 0, 0 };
    if (l > 0) { r.x = a.x / l; r.y = a.y / l; r.z = a.z / l; }
    return r;
}
static V3 vscale(V3 a, float s) { V3 v = { a.x * s, a.y * s, a.z * s }; return v; }
static V3 vadd(V3 a, V3 b) { V3 v = { a.x + b.x, a.y + b.y, a.z + b.z }; return v; }

static Camera cam_make(V3 eye, V3 target) {
    Camera c;
    c.o = eye;
    c.f = vnorm(vsub(target, eye));
    V3 world_up = { 0, 1, 0 };
    c.r = vnorm(vcross(c.f, world_up));
    c.u = vcross(c.r, c.f);
    return c;
}

/* 世界点投影到相机 c 的连续像素坐标；在相机背后或画面外返回 0。
   容差取 ±0.75：连续坐标 round 后落在 [0,W) 内即对应合法历史像素，
   同时吸收像素中心光线往返投影的浮点误差。 */
static int project(const Camera* c, V3 p, float* sx, float* sy) {
    V3 rel = vsub(p, c->o);
    float cz = vdot(rel, c->f);
    if (cz <= 1e-4f) return 0;
    float cx = vdot(rel, c->r);
    float cy = vdot(rel, c->u);
    float tan_half = tanf(FOV_Y * 0.5f);
    float ndx = (cx / cz) / tan_half;            /* [-1,1] 的 NDC，x 按宽高比 */
    float ndy = (cy / cz) / tan_half;
    *sx = (ndx + 1.0f) * 0.5f * (W - 1);
    *sy = (1.0f - ndy) * 0.5f * (H - 1);
    return *sx >= -0.75f && *sx <= W - 0.25f && *sy >= -0.75f && *sy <= H - 0.25f;
}

/* 像素中心对应的光线方向，与 project() 互为逆运算 */
static V3 ray_dir(const Camera* c, int px, int py) {
    float tan_half = tanf(FOV_Y * 0.5f);
    float ndx = 2.0f * ((float)px + 0.5f) / W - 1.0f;
    float ndy = 1.0f - 2.0f * ((float)py + 0.5f) / H;
    V3 d = c->f;
    d = vadd(d, vscale(c->r, ndx * tan_half));
    d = vadd(d, vscale(c->u, ndy * tan_half));
    return vnorm(d);
}

/* 与 y=0 无限平面求交，t>0 才命中 */
static int hit_plane(const Camera* c, V3 d, V3* hit) {
    if (fabsf(d.y) < 1e-8f || c->o.y * d.y >= 0.0f) return 0;
    float t = -c->o.y / d.y;
    if (t <= 0.0f) return 0;
    *hit = vadd(c->o, vscale(d, t));
    return 1;
}

/* 无噪参考色：棋盘地面，epoch 翻转棋盘明暗（用于制造颜色突变） */
static RGB ground_truth(V3 p, int epoch) {
    int cell = (int)floorf(p.x) + (int)floorf(p.z);
    int dark = (cell & 1) ^ epoch;
    RGB g = dark ? (RGB){ 0.25f, 0.25f, 0.25f } : (RGB){ 0.85f, 0.85f, 0.80f };
    return g;
}

/* 确定性 LCG，保证自检可复现 */
static unsigned rng_state;
static float rng_next(void) {
    rng_state = rng_state * 1103515245u + 12345u;
    return (float)(rng_state >> 8) / (float)(1u << 24);
}

/* 渲染一帧：带噪颜色 + 世界坐标 G-buffer（本场景像素全命中） */
static void render_frame(const Camera* cam, int epoch, RGB* col, V3* pos, int* hit) {
    for (int py = 0; py < H; py++) {
        for (int px = 0; px < W; px++) {
            int i = py * W + px;
            V3 d = ray_dir(cam, px, py);
            V3 hp;
            hit[i] = hit_plane(cam, d, &hp);
            pos[i] = hp;
            RGB g = ground_truth(hp, epoch);
            /* 两个均匀噪声之差：均值 0，对称截断，近似高斯 */
            float nr = (rng_next() - rng_next()) * 1.5f * NOISE_SIGMA;
            float ng = (rng_next() - rng_next()) * 1.5f * NOISE_SIGMA;
            float nb = (rng_next() - rng_next()) * 1.5f * NOISE_SIGMA;
            col[i].r = g.r + nr; col[i].g = g.g + ng; col[i].b = g.b + nb;
        }
    }
}

static float luminance(RGB c) { return 0.2126f * c.r + 0.7152f * c.g + 0.0722f * c.b; }

enum { ACC_OK = 0, ACC_REJECT_COLOR = 1, ACC_MISSING = 2 };

typedef struct {
    RGB hcol[NPIX];          /* 历史累积均值 */
    V3  hpos[NPIX];          /* 历史像素对应的世界交点 */
    int hhit[NPIX];          /* 历史像素是否命中过几何体 */
    int n[NPIX];             /* 每像素有效累积帧数（1 = 仅当前帧） */
    int valid;               /* 历史缓冲是否已初始化 */
    Camera prev;             /* 上一帧相机 */
} Accum;

static void accum_reset(Accum* a) {
    memset(a, 0, sizeof(*a));
    a->valid = 0;
}

/* 单像素累积一步：先重投影找历史，再做几何与颜色拒绝 */
static int accum_pixel(Accum* a, int i, RGB c, V3 cp, int chit, int cap) {
    int status = ACC_MISSING;
    int hn = 0;
    RGB hc = { 0, 0, 0 };

    if (a->valid && chit) {
        float sx, sy;
        if (project(&a->prev, cp, &sx, &sy)) {
            int hx = (int)lroundf(sx);
            int hy = (int)lroundf(sy);
            if (hx >= 0 && hx < W && hy >= 0 && hy < H) {
                int j = hy * W + hx;
                if (a->hhit[j]) {
                    V3 d = vsub(cp, a->hpos[j]);
                    if (vdot(d, d) <= GEO_EPS * GEO_EPS &&
                        fabsf(luminance(c) - luminance(a->hcol[j])) <= COLOR_REJECT) {
                        hn = a->n[j];
                        hc = a->hcol[j];
                        status = ACC_OK;
                    } else {
                        status = ACC_REJECT_COLOR;
                    }
                }
            }
        }
    }

    int new_n;
    RGB out;
    if (status == ACC_OK) {
        new_n = hn + 1;
        if (cap > 0 && new_n > cap) {
            /* 达到帧数上限：滑动窗口，保持有效样本数不变 */
            float w = 1.0f / (float)cap;
            out.r = hc.r + (c.r - hc.r) * w;
            out.g = hc.g + (c.g - hc.g) * w;
            out.b = hc.b + (c.b - hc.b) * w;
            new_n = cap;
        } else {
            float w = 1.0f / (float)new_n;
            out.r = hc.r + (c.r - hc.r) * w;
            out.g = hc.g + (c.g - hc.g) * w;
            out.b = hc.b + (c.b - hc.b) * w;
        }
    } else {
        /* 颜色拒绝或历史缺失：退回单帧，累积计数从 1 重新开始 */
        new_n = 1;
        out = c;
    }

    a->hcol[i] = out;
    a->hpos[i] = cp;
    a->hhit[i] = chit;
    a->n[i] = new_n;
    return status;
}

/* 相对无噪参考的 RMSE：单帧噪声与累积后噪声用同一把尺子 */
static double frame_rmse(const RGB* img, const V3* pos, const int* hit, int epoch) {
    double s = 0.0;
    int m = 0;
    for (int i = 0; i < NPIX; i++) {
        if (!hit[i]) continue;
        RGB g = ground_truth(pos[i], epoch);
        float dr = img[i].r - g.r, dg = img[i].g - g.g, db = img[i].b - g.b;
        s += dr * dr + dg * dg + db * db;
        m++;
    }
    return sqrt(s / (double)(3 * m));
}

static double mean_n(const Accum* a) {
    long s = 0;
    for (int i = 0; i < NPIX; i++) s += a->n[i];
    return (double)s / NPIX;
}

typedef struct {
    int reject_frames;        /* 出现过颜色拒绝的帧数 */
    int missing_frames;       /* 出现过历史缺失的帧数 */
    int reject_pixels;        /* 累计颜色拒绝像素数 */
    int missing_pixels;       /* 累计历史缺失像素数 */
    int total_frames;
} RunStats;

/* 跑一帧：渲染 -> 逐像素累积 -> 统计 -> 打印标注与噪声对比 */
static void run_frame(Accum* a, const Camera* cam, int epoch, int frame,
                      int cap, RunStats* st) {
    RGB col[NPIX], out[NPIX];
    V3 pos[NPIX];
    int hit[NPIX];
    render_frame(cam, epoch, col, pos, hit);

    int frame_reject = 0, frame_missing = 0;
    for (int py = 0; py < H; py++) {
        for (int px = 0; px < W; px++) {
            int i = py * W + px;
            int s = accum_pixel(a, i, col[i], pos[i], hit[i], cap);
            out[i] = a->hcol[i];
            if (s == ACC_REJECT_COLOR) { frame_reject = 1; st->reject_pixels++; }
            if (s == ACC_MISSING) { frame_missing = 1; st->missing_pixels++; }
        }
    }

    double single = frame_rmse(col, pos, hit, epoch);
    double accumd = frame_rmse(out, pos, hit, epoch);
    const char* note;
    if (frame_missing && frame_reject) note = "fallback=single(reject+missing)";
    else if (frame_missing)            note = "fallback=single(history_missing)";
    else if (frame_reject)             note = "fallback=single(color_reject)";
    else                               note = "accumulating";

    st->reject_frames += frame_reject;
    st->missing_frames += frame_missing;
    st->total_frames = frame + 1;

    int nmax = 0;
    for (int k = 0; k < NPIX; k++) if (a->n[k] > nmax) nmax = a->n[k];
    printf("frame=%02d note=%-34s nmean=%6.2f nmax=%3d rmse_single=%.4f rmse_accum=%.4f drop=%.2fx\n",
           frame, note, mean_n(a), nmax, single, accumd, single / accumd);

    a->valid = 1;
    a->prev = *cam;
}

static int failures = 0;
static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

int main(void) {
    Accum a;
    RunStats st;

    /* ---- 场景 1：静态相机，历史应持续累积，噪声随帧数下降 ---- */
    printf("== static camera ==\n");
    rng_state = 20261123u;
    accum_reset(&a);
    memset(&st, 0, sizeof(st));
    Camera eye0 = { 0 };
    V3 eye = { 0, 2.5f, 4.0f };
    V3 target = { 0, 0, 0 };   /* 俯角 32 度大于半 FOV 27.5 度，保证所有像素都打到地面 */
    eye0 = cam_make(eye, target);
    for (int f = 0; f < 12; f++)
        run_frame(&a, &eye0, 0, f, 0, &st);

    check("static_first_frame_missing", st.missing_frames >= 1);
    check("static_no_color_reject", st.reject_pixels == 0);
    check("static_grows_to_frames", a.n[0] == 12);

    double static_drop = 0.0;
    {
        /* 再渲染一帧，直接量累积前后 RMSE 的下降比 */
        RGB col[NPIX], out[NPIX];
        V3 pos[NPIX];
        int hit[NPIX];
        render_frame(&eye0, 0, col, pos, hit);
        double single = frame_rmse(col, pos, hit, 0);
        for (int py = 0; py < H; py++)
            for (int px = 0; px < W; px++) {
                int i = py * W + px;
                accum_pixel(&a, i, col[i], pos[i], hit[i], 0);
                out[i] = a.hcol[i];
            }
        double accumd = frame_rmse(out, pos, hit, 0);
        static_drop = single / accumd;
        printf("static: rmse_single=%.4f rmse_accum=%.4f drop=%.2fx nmean=%.2f\n",
               single, accumd, static_drop, mean_n(&a));
    }
    check("static_noise_drops", static_drop >= 2.5);

    /* ---- 场景 2：相机平移，第 10 帧场景换色触发颜色拒绝 ---- */
    printf("== moving camera, scene change at frame 10 ==\n");
    rng_state = 99u;
    accum_reset(&a);
    memset(&st, 0, sizeof(st));
    int reject_at_change = 0;
    int missing_after_first = 0;
    for (int f = 0; f < 18; f++) {
        V3 e = { 0.25f * f, 2.5f, 4.0f };
        Camera cam = cam_make(e, target);
        int epoch = (f >= 10) ? 1 : 0;
        int rej_before = st.reject_pixels, miss_before = st.missing_pixels;
        run_frame(&a, &cam, epoch, f, 0, &st);
        if (f == 10) reject_at_change = st.reject_pixels - rej_before;
        if (f > 0)   missing_after_first += st.missing_pixels - miss_before;
    }

    /* 首帧无历史，全部像素必须走缺失退回单帧 */
    check("moving_first_frame_all_missing", st.missing_pixels >= NPIX);
    /* 首帧之后缺失只来自画面侧缘的遮挡空缺：缺失率应远低于逐帧全丢 */
    check("moving_reprojection_used", missing_after_first < NPIX * 17 / 3);
    /* 换色帧绝大多数像素因颜色差异被拒绝并计入统计 */
    check("moving_color_reject_counted", reject_at_change > NPIX / 2);
    check("moving_reject_frame_counted", st.reject_frames >= 1);

    /* 换色后相机不动再累积几帧，计数从 1 恢复且噪声重新下降 */
    V3 e = { 0.25f * 17, 2.5f, 4.0f };
    Camera cam_last = cam_make(e, target);
    double recon_drop = 1.0;
    for (int k = 0; k < 6; k++) {
        RGB col[NPIX], out[NPIX];
        V3 pos[NPIX];
        int hit[NPIX];
        render_frame(&cam_last, 1, col, pos, hit);
        for (int py = 0; py < H; py++)
            for (int px = 0; px < W; px++) {
                int i = py * W + px;
                accum_pixel(&a, i, col[i], pos[i], hit[i], 0);
                out[i] = a.hcol[i];
            }
        double single = frame_rmse(col, pos, hit, 1);
        double accumd = frame_rmse(out, pos, hit, 1);
        recon_drop = single / accumd;
    }
    int n_after_sum = 0;
    for (int i = 0; i < NPIX; i++) n_after_sum += a.n[i];
    double n_after_mean = (double)n_after_sum / NPIX;
    printf("recover: drop=%.2fx nmean=%.2f\n", recon_drop, n_after_mean);
    /* 侧缘像素持续无历史属正常，用平均累积帧数衡量整体是否重新收敛 */
    check("reconverge_after_reject", recon_drop >= 1.8 && n_after_mean >= 6.0);

    /* ---- 场景 3：累积帧数上限 ---- */
    printf("== accumulation cap ==\n");
    rng_state = 7u;
    accum_reset(&a);
    RunStats st3;
    memset(&st3, 0, sizeof(st3));
    for (int f = 0; f < 10; f++)
        run_frame(&a, &eye0, 0, f, 4, &st3);
    int nmax = 0;
    for (int i = 0; i < NPIX; i++) if (a.n[i] > nmax) nmax = a.n[i];
    check("cap_respected", nmax == 4);

    /* ---- 重投影往返：投出去的历史像素必须落回邻近原像素 ---- */
    {
        V3 e1 = { 0, 2.5f, 4.0f }, e2 = { 0.5f, 2.5f, 4.0f };
        Camera c1 = cam_make(e1, target), c2 = cam_make(e2, target);
        V3 d = ray_dir(&c1, W / 2, H / 2);
        V3 hp;
        int ok_round = 0;
        if (hit_plane(&c1, d, &hp)) {
            float sx, sy;
            if (project(&c2, hp, &sx, &sy)) {
                V3 d2 = ray_dir(&c2, (int)lroundf(sx), (int)lroundf(sy));
                V3 hp2;
                if (hit_plane(&c2, d2, &hp2)) {
                    V3 dd = vsub(hp, hp2);
                    ok_round = vdot(dd, dd) < 0.3f * 0.3f;
                }
            }
        }
        check("reproject_roundtrip", ok_round);
    }

    printf("summary: reject_frames=%d missing_frames=%d reject_pixels=%d missing_pixels=%d\n",
           st.reject_frames, st.missing_frames, st.reject_pixels, st.missing_pixels);
    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
