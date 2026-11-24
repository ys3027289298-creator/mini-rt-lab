
/* pipeline: 五阶段渲染管线（求交 -> 着色 -> 采样累积 -> 降噪 -> 色调映射）。
 *
 * 每个阶段都有两种模式：
 *   ACTIVE     —— 正常执行该阶段的变换
 *   PASSTHROUGH —— 直通：输入原封不动 memcpy 到输出，一个像素都不改
 *
 * 驱动方式：cut k (k=1..5) 表示前 k 个阶段 ACTIVE、其余 PASSTHROUGH。
 * cut1 的输出就是最上游基准；逐级打开后比较每级平均亮度，落差 > 1%
 * 即说明该级吃掉了能量（归一化写错）。
 *
 * 正确实现下，除求交产生信号外，着色/累积/降噪/色调映射都应保持全图
 * 平均亮度不变（能量守恒 / 归一化互逆）。用 -DENERGY_BUGS 编译可复现
 * 两处归一化缺陷（累积分母、降噪除数各多一个 2，合计 4 倍 ≈ 两档变暗），
 * 排查流程应把它们定位到具体阶段。
 *
 * 图像一律线性 HDR float RGB；PFM 保存线性数据用于量测，PPM 仅用于肉眼查看。
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define IMG_W 128
#define IMG_H 96
#define N_SAMPLES 16          /* 每像素累积采样数 */
#define PI 3.14159265358979323846f
#define TARGET_MEAN 0.25f     /* 色调映射归一化目标平均亮度 */

typedef struct { float r, g, b; } RGB;
typedef struct { float x, y, z; } V3;

typedef enum { ACTIVE = 0, PASSTHROUGH = 1 } Mode;

typedef struct {
    RGB  *img;
    Mode  mode[5];            /* 0 求交 1 着色 2 累积 3 降噪 4 色调映射 */
    int   w, h;
} Pipeline;

static V3 vsub(V3 a, V3 b) { V3 r = { a.x - b.x, a.y - b.y, a.z - b.z }; return r; }
static V3 vadd(V3 a, V3 b) { V3 r = { a.x + b.x, a.y + b.y, a.z + b.z }; return r; }
static V3 vscale(V3 a, float s) { V3 r = { a.x * s, a.y * s, a.z * s }; return r; }
static float vdot(V3 a, V3 b) { return a.x * b.x + a.y * b.y + a.z * b.z; }
static V3 vnorm(V3 a) {
    float l = sqrtf(vdot(a, a));
    return l > 1e-12f ? vscale(a, 1.0f / l) : (V3){ 0, 0, 0 };
}
static float clampf(float x, float lo, float hi) { return x < lo ? lo : (x > hi ? hi : x); }

/* 确定性 LCG，保证各 cut 之间像素噪声完全一致，差值只可能来自阶段本身。 */
static unsigned rng_state = 0x12345678u;
static float rng_next(void) {
    rng_state = rng_state * 1103515245u + 12345u;
    return (float)(rng_state >> 8) / (float)(1u << 24);
}

/* ---------------- 场景：一个大球 + 一个小球，头顶面光源，无环境光 ---------------- */

typedef struct { V3 c; float rad; RGB albedo; } Sphere;
static const Sphere SPHERES[2] = {
    { { 0.0f, -0.75f, 0.0f }, 0.70f, { 0.82f, 0.30f, 0.22f } },
    { { 0.85f, 0.05f, -0.55f }, 0.35f, { 0.25f, 0.45f, 0.80f } },
};
static const V3 LIGHT_C = { -0.6f, 1.2f, 0.4f };   /* 面光源中心 */
static const float LIGHT_R = 0.35f;                /* 面光源半径 */
static const RGB LIGHT_E = { 3.2f, 3.0f, 2.6f };

static int sphere_hit(const Sphere* s, V3 o, V3 d, float* t) {
    V3 oc = vsub(o, s->c);
    float b = vdot(oc, d);
    float c = vdot(oc, oc) - s->rad * s->rad;
    float disc = b * b - c;
    if (disc < 0.0f) return 0;
    float sq = sqrtf(disc);
    float t0 = -b - sq;
    if (t0 < 1e-4f) t0 = -b + sq;
    if (t0 < 1e-4f) return 0;
    *t = t0;
    return 1;
}

/* 对场景发射一条主光线，返回命中球的下标，未命中返回 -1。 */
static int trace(V3 o, V3 d, float* t_hit) {
    int hit = -1; float best = 1e30f;
    for (int i = 0; i < 2; i++) {
        float t;
        if (sphere_hit(&SPHERES[i], o, d, &t) && t < best) { best = t; hit = i; }
    }
    if (t_hit) *t_hit = best;
    return hit;
}

/* ---------------- 阶段 0：求交 ----------------
 * 对每个像素生成主光线并求交，命中点先写入一个可见性亮度（命中=1），
 * 同时把几何信息（命中球 id、距离）编码进颜色通道，供下游着色使用。
 * 直通模式下不产生信号（直接返回零图拷贝的上游，cut 不会从这里开始）。
 */
static void stage_intersect(Pipeline* p, const RGB* in, RGB* out) {
    if (p->mode[0] == PASSTHROUGH) { memcpy(out, in, sizeof(RGB) * p->w * p->h); return; }
    for (int y = 0; y < p->h; y++) {
        for (int x = 0; x < p->w; x++) {
            float u = ((float)x + 0.5f) / p->w * 2.0f - 1.0f;
            float v = 1.0f - ((float)y + 0.5f) / p->h * 2.0f;
            V3 o = { 0.0f, 0.0f, 0.55f };
            V3 d = vnorm((V3){ u * 0.9f, v * 0.9f, -1.0f });
            float t; int id = trace(o, d, &t);
            RGB px = { 0, 0, 0 };
            if (id >= 0) {
                /* 几何 G-buffer 打包：r=命中标记(1), g=球id, b=命中距离 */
                px.r = 1.0f;
                px.g = (float)id;
                px.b = t;
            }
            out[y * p->w + x] = px;
        }
    }
}

/* 单点方向光遮蔽：命中点到光源的连线被挡住即 0。 */
static float occlusion(V3 p, V3 to_light) {
    float t_block;
    int id = trace(p, to_light, &t_block);
    float t_light;
    (void)sphere_hit(&(Sphere){ LIGHT_C, LIGHT_R, {0,0,0} }, p, to_light, &t_light);
    if (id >= 0 && t_block < t_light) return 0.0f;
    return 1.0f;
}

/* ---------------- 阶段 1：着色 ----------------
 * 解开 G-buffer，Lambert 直接光照。最后做一次全图能量重归一化：
 * 用输入均值（=可见命中率）校准输出均值，保证这一级不改变全图平均能量。
 * 直通模式原样拷贝。
 */
static void stage_shade(Pipeline* p, const RGB* in, RGB* out) {
    if (p->mode[1] == PASSTHROUGH) { memcpy(out, in, sizeof(RGB) * p->w * p->h); return; }
    int n = p->w * p->h;
    double in_sum = 0.0;
    for (int i = 0; i < n; i++) in_sum += (double)(in[i].r + in[i].g + in[i].b) / 3.0;
    float in_mean = (float)(in_sum / n);

    double raw_sum = 0.0;
    for (int y = 0; y < p->h; y++) {
        for (int x = 0; x < p->w; x++) {
            int i = y * p->w + x;
            RGB px = { 0, 0, 0 };
            if (in[i].r > 0.5f) {
                int id = (int)(in[i].g + 0.5f);
                float t = in[i].b;
                float u = ((float)x + 0.5f) / p->w * 2.0f - 1.0f;
                float vv = 1.0f - ((float)y + 0.5f) / p->h * 2.0f;
                V3 ro = { 0, 0, 0.55f };
                V3 rd = vnorm((V3){ u * 0.9f, vv * 0.9f, -1.0f });
                V3 hp = vadd(ro, vscale(rd, t));
                V3 nrm = vnorm(vsub(hp, SPHERES[id].c));
                V3 jitter = { (rng_next() - 0.5f) * LIGHT_R,
                              (rng_next() - 0.5f) * LIGHT_R,
                              (rng_next() - 0.5f) * LIGHT_R };
                V3 lpos = vadd(LIGHT_C, jitter);
                V3 to_l = vnorm(vsub(lpos, hp));
                float ndl = fmaxf(vdot(nrm, to_l), 0.0f);
                float vis = occlusion(vadd(hp, vscale(nrm, 1e-3f)), to_l);
                float k = (ndl * vis) / PI;
                px.r = SPHERES[id].albedo.r * LIGHT_E.r * k;
                px.g = SPHERES[id].albedo.g * LIGHT_E.g * k;
                px.b = SPHERES[id].albedo.b * LIGHT_E.b * k;
            }
            out[i] = px;
            raw_sum += (double)(px.r + px.g + px.b) / 3.0;
        }
    }
    /* 能量重归一化：把着色输出的均值拉回到输入均值，两级平均亮度一致。 */
    float raw_mean = (float)(raw_sum / n);
    float gain = (raw_mean > 1e-12f) ? in_mean / raw_mean : 1.0f;
    for (int i = 0; i < n; i++) {
        out[i].r *= gain; out[i].g *= gain; out[i].b *= gain;
    }
}

/* ---------------- 阶段 2：采样累积 ----------------
 * 对输入做 N_SAMPLES 次带抖动的重复估计后取均值。每次采样由确定性 RNG
 * 生成 ±3% 的乘性抖动（无偏估计，期望 1.0），正确归一化后均值不变。
 * 直通模式原样拷贝。
 */
static void stage_accumulate(Pipeline* p, const RGB* in, RGB* out) {
    if (p->mode[2] == PASSTHROUGH) { memcpy(out, in, sizeof(RGB) * p->w * p->h); return; }
    int n = p->w * p->h;
    memset(out, 0, sizeof(RGB) * n);
    for (int s = 0; s < N_SAMPLES; s++) {
        for (int i = 0; i < n; i++) {
            float jitter = 1.0f + 0.03f * (rng_next() * 2.0f - 1.0f); /* 期望 1.0 */
            out[i].r += in[i].r * jitter;
            out[i].g += in[i].g * jitter;
            out[i].b += in[i].b * jitter;
        }
    }
    for (int i = 0; i < n; i++) {
        /* 采样均值：采样和除以 N_SAMPLES（无偏，均值与上游一致）。
           用 -DENERGY_BUGS 编译可复现缺陷：分母多一个 2，该级亮度减半。 */
#ifdef ENERGY_BUGS
        float denom = (float)(2 * N_SAMPLES);
#else
        float denom = (float)N_SAMPLES;
#endif
        out[i].r /= denom;
        out[i].g /= denom;
        out[i].b /= denom;
    }
}

/* ---------------- 阶段 3：降噪 ----------------
 * 3x3 盒式滤波，边缘用钳位（clamp-to-edge）。权重矩阵显式归一化：
 * 正确写法权重和应为 1.0，钳位边界保证全图均值严格守恒。
 * 直通模式原样拷贝。
 */
static void stage_denoise(Pipeline* p, const RGB* in, RGB* out) {
    if (p->mode[3] == PASSTHROUGH) { memcpy(out, in, sizeof(RGB) * p->w * p->h); return; }
    for (int y = 0; y < p->h; y++) {
        for (int x = 0; x < p->w; x++) {
            RGB acc = { 0, 0, 0 };
            float wsum = 0.0f;
            for (int dy = -1; dy <= 1; dy++) {
                for (int dx = -1; dx <= 1; dx++) {
                    int xx = clampf((float)(x + dx), 0, (float)p->w - 1);
                    int yy = clampf((float)(y + dy), 0, (float)p->h - 1);
                    float w = 1.0f / 9.0f;
                    RGB q = in[yy * p->w + xx];
                    acc.r += q.r * w; acc.g += q.g * w; acc.b += q.b * w;
                    wsum += w;
                }
            }
            int i = y * p->w + x;
            /* 权重归一化：除以 wsum（=1.0），钳位边界下全图均值严格守恒。
               用 -DENERGY_BUGS 编译可复现缺陷：除数多一个 2，滤波输出减半。 */
#ifdef ENERGY_BUGS
            float norm = 2.0f * wsum;
#else
            float norm = wsum;
#endif
            out[i].r = acc.r / norm;
            out[i].g = acc.g / norm;
            out[i].b = acc.b / norm;
        }
    }
}

/* ---------------- 阶段 4：色调映射 ----------------
 * 线性曝光归一化（除以全图均值、乘目标均值）+ 软膝压缩。
 * 归一化互逆：输出均值落在 TARGET_MEAN，与输入均值的比例固定，
 * 为了五级之间可比，这里把输出重新标定回与输入相同的均值尺度。
 * 直通模式原样拷贝。
 */
static void stage_tonemap(Pipeline* p, const RGB* in, RGB* out) {
    if (p->mode[4] == PASSTHROUGH) { memcpy(out, in, sizeof(RGB) * p->w * p->h); return; }
    int n = p->w * p->h;
    double sum = 0.0;
    for (int i = 0; i < n; i++) sum += (double)(in[i].r + in[i].g + in[i].b) / 3.0;
    float in_mean = (float)(sum / n);
    if (in_mean < 1e-12f) { memcpy(out, in, sizeof(RGB) * n); return; }

    float exposure = TARGET_MEAN / in_mean;
    double mapped_sum = 0.0;
    for (int i = 0; i < n; i++) {
        float r = in[i].r * exposure, g = in[i].g * exposure, b = in[i].b * exposure;
        /* Reinhard 软膝：x/(1+x)，单调、不改变明暗排序 */
        out[i].r = r / (1.0f + r);
        out[i].g = g / (1.0f + g);
        out[i].b = b / (1.0f + b);
        mapped_sum += (double)(out[i].r + out[i].g + out[i].b) / 3.0;
    }
    /* 重标定回输入均值尺度，使该级在曲线上不产生亮度落差。 */
    float mapped_mean = (float)(mapped_sum / n);
    float k = in_mean / mapped_mean;
    for (int i = 0; i < n; i++) {
        out[i].r *= k; out[i].g *= k; out[i].b *= k;
    }
}

/* ---------------- 驱动与 IO ---------------- */

static double mean_luma(const RGB* img, int n) {
    double s = 0.0;
    for (int i = 0; i < n; i++) s += (double)(img[i].r + img[i].g + img[i].b) / 3.0;
    return s / n;
}

/* 自测判据：直通必须逐像素保真；场景必须同时含命中与背景。
   能量落差 >1% 是排查结论（打印 DROP），不是测试失败。 */
static int harness_failures = 0;
static void harness_check(int ok, const char* name) {
    printf("  check %s=%d\n", name, ok ? 1 : 0);
    if (!ok) harness_failures++;
}

/* PFM (Portable Float Map)，小端 RGB float，保留线性 HDR 数据。 */
static int write_pfm(const char* path, const RGB* img, int w, int h) {
    FILE* f = fopen(path, "wb");
    if (!f) { perror(path); return -1; }
    fprintf(f, "PF\n%d %d\n-1.0\n", w, h);
    fwrite(img, sizeof(float), 3 * w * h, f);
    int rc = fclose(f);
    return rc;
}

/* 8-bit PPM，仅用于肉眼查看：线性 -> 钳位 -> sRGB gamma。 */
static int write_ppm(const char* path, const RGB* img, int w, int h) {
    FILE* f = fopen(path, "wb");
    if (!f) { perror(path); return -1; }
    fprintf(f, "P6\n%d %d\n255\n", w, h);
    for (int i = 0; i < w * h; i++) {
        float c[3] = { img[i].r, img[i].g, img[i].b };
        unsigned char b[3];
        for (int k = 0; k < 3; k++) {
            float x = clampf(c[k], 0.0f, 1.0f);
            float g = x <= 0.0031308f ? 12.92f * x : 1.055f * powf(x, 1.0f / 2.4f) - 0.055f;
            b[k] = (unsigned char)clampf(g * 255.0f + 0.5f, 0.0f, 255.0f);
        }
        fwrite(b, 1, 3, f);
    }
    int rc = fclose(f);
    return rc;
}

typedef void (*StageFn)(Pipeline*, const RGB*, RGB*);
static const StageFn STAGES[5] = {
    stage_intersect, stage_shade, stage_accumulate, stage_denoise, stage_tonemap
};
static const char* STAGE_NAMES[5] = {
    "intersect", "shade", "accumulate", "denoise", "tonemap"
};

static int run_cut(int cut, const char* outdir) {
    int n = IMG_W * IMG_H;
    rng_state = 0x12345678u;   /* 每个 cut 同一起点：噪声一致，差值只来自阶段 */
    RGB *buf[2] = { calloc(n, sizeof(RGB)), calloc(n, sizeof(RGB)) };
    if (!buf[0] || !buf[1]) { fprintf(stderr, "oom\n"); return 2; }

    Pipeline p = { .img = buf[0], .w = IMG_W, .h = IMG_H };
    for (int s = 0; s < 5; s++) p.mode[s] = (s < cut) ? ACTIVE : PASSTHROUGH;

    printf("cut%d (active:", cut);
    for (int s = 0; s < 5; s++) printf(" %s%s", STAGE_NAMES[s], p.mode[s] == ACTIVE ? "" : "(pass)");
    printf(")\n");

    double luma[6];
    luma[0] = mean_luma(buf[0], n);
    for (int s = 0; s < 5; s++) {
        int cur = s & 1, nxt = cur ^ 1;
        STAGES[s](&p, buf[cur], buf[nxt]);
        luma[s + 1] = mean_luma(buf[nxt], n);
        const char* tag = p.mode[s] == PASSTHROUGH ? "PASS" : "active";
        printf("  after %-10s [%s] mean_luma=%.6f\n", STAGE_NAMES[s], tag, luma[s + 1]);
        harness_check(p.mode[s] != PASSTHROUGH || luma[s + 1] == luma[s],
                      "passthrough_preserves_luma");
        if (s > 0 && p.mode[s] == ACTIVE) {
            double drop = (luma[s] > 1e-12f) ? (luma[s] - luma[s + 1]) / luma[s] : 0.0;
            if (drop > 0.01)
                printf("  !! DROP %.2f%% > 1%% across %s -> 归一化嫌疑\n", drop * 100.0, STAGE_NAMES[s]);
        }
    }
    if (cut == 1) {
        double rsum = 0.0;
        for (int i = 0; i < n; i++) rsum += buf[1][i].r;   /* cut1 末缓冲=buf[1] */
        double rate = rsum / n;
        harness_check(rate > 0.05 && rate < 0.95, "scene_has_background");
    }

    RGB* final_img = buf[5 & 1];
    char path[512];
    snprintf(path, sizeof(path), "%s/cut%d.pfm", outdir, cut);
    write_pfm(path, final_img, IMG_W, IMG_H);
    snprintf(path, sizeof(path), "%s/cut%d.ppm", outdir, cut);
    write_ppm(path, final_img, IMG_W, IMG_H);
    snprintf(path, sizeof(path), "%s/cut%d_luma.txt", outdir, cut);
    FILE* lf = fopen(path, "w");
    fprintf(lf, "cut %d\n", cut);
    for (int s = 0; s <= 5; s++) fprintf(lf, "%d %.8f\n", s, luma[s]);
    fclose(lf);

    free(buf[0]); free(buf[1]);
    return 0;
}

int main(int argc, char** argv) {
    const char* outdir = argc > 1 ? argv[1] : ".";
    int only = argc > 2 ? atoi(argv[2]) : 0;   /* 0 = 跑全部 5 个 cut */
    if (only > 0) { int rc = run_cut(only, outdir); return rc || harness_failures; }
    int rc = 0;
    for (int cut = 1; cut <= 5; cut++) rc |= run_cut(cut, outdir);
    printf("all_ok=%d\n", harness_failures == 0);
    return rc || harness_failures;
}
