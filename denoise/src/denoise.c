
/* denoise: 采样器/降噪器统一接口与时间累积。
 *
 * 复现“接口对不上 -> 画面忽明忽暗”：
 *   1) 同一个确定性噪声序列（xorshift64 + Box-Muller，固定种子）同时喂给
 *      旧路径和统一路径，保证两路径输入逐样本一致；
 *   2) 旧路径采样器每帧 spp 在 4/8 间交替，降噪器按自己缓存的帧计数做固定
 *      权重时间混合，并把“均值”再除一次 spp、把单样本方差当均值方差用——
 *      spp 错位让输出随帧 2:1 交替变亮变暗；
 *   3) 统一路径走 rt_frame.h 契约，按总样本数加权，闪烁帧数必须为 0。
 *
 * 自检输出 check 风格与 bvh/sample/material 三个模块一致。
 */
#include "rt_frame.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>

#define W 64
#define H 48
#define PX (W * H)
#define FRAMES 40
#define WARMUP 12                 /* 统计闪烁前丢弃的暖机帧 */
#define FLICKER_EPS 2e-3          /* 帧间均值相对变化阈值 */
#define SEED 0x243f6a8885a308d3ULL
#define GOLDEN_MAGIC 0x52544731u /* "RTG1" */

/* ---------------- RNG 与确定性场景 ---------------- */

static unsigned long long rng_state;
static unsigned long long xs64(void) {
    unsigned long long x = rng_state;
    x ^= x << 13; x ^= x >> 7; x ^= x << 17;
    rng_state = x;
    return x;
}
static double rand_u(void) { return (double)(xs64() >> 11) / (double)(1ULL << 53); }

static double gauss(void) {
    double u1 = rand_u(); if (u1 < 1e-12) u1 = 1e-12;
    double u2 = rand_u();
    return sqrt(-2.0 * log(u1)) * cos(2.0 * 3.14159265358979323846 * u2);
}

/* 确定性“地面真值”：两个高斯光斑 + 斜坡 + 暗角，全图均值约 0.25 */
static double ground_truth(int x, int y) {
    double nx = (double)x / (W - 1), ny = (double)y / (H - 1);
    double g1 = 0.55 * exp(-(((nx - 0.30) * (nx - 0.30)) + ((ny - 0.35) * (ny - 0.35))) / 0.020);
    double g2 = 0.35 * exp(-(((nx - 0.72) * (nx - 0.72)) + ((ny - 0.65) * (ny - 0.65))) / 0.035);
    double ramp = 0.12 * (nx + ny) * 0.5;
    double vig = 0.9 + 0.1 * cos(3.14159265 * (nx - 0.5));
    double v = (0.10 + g1 + g2 + ramp) * vig;
    return v < 0 ? 0 : v;
}

/* 单像素单样本：真值 + 幅度随真值变化的高斯噪声（暗区噪点更密，符合 MC 特性） */
static float take_sample(int x, int y) {
    double gt = ground_truth(x, y);
    double sigma = 0.30 * sqrt(gt + 0.02);
    double v = gt + sigma * gauss();
    return (float)(v < 0 ? 0 : v);
}

static int spp_of(int frame) { return (frame & 1) ? 4 : 8; } /* 4/8 交替 */

/* ---------------- 统一接口实现 ---------------- */

int rt_accum_init(RtAccum* a, int width, int height) {
    memset(a, 0, sizeof(*a));
    a->width = width; a->height = height;
    size_t n = (size_t)width * height;
    a->mean = malloc(n * sizeof(float));
    a->vpool = calloc(n, sizeof(float));
    a->mvar = calloc(n, sizeof(float));
    if (!a->mean || !a->vpool || !a->mvar) {
        rt_accum_free(a);
        return -1;
    }
    return 0;
}

void rt_accum_free(RtAccum* a) {
    free(a->mean); free(a->vpool); free(a->mvar);
    memset(a, 0, sizeof(*a));
}

int rt_accum_push(RtAccum* a, const RtFrame* f) {
    if (f->width != a->width || f->height != a->height || f->spp <= 0) return -1;
    long N_old = a->samples;
    long N_new = N_old + f->spp;
    double inv = 1.0 / (double)N_new;
    size_t npix = (size_t)a->width * a->height;
    for (size_t i = 0; i < npix; i++) {
        /* 按样本数加权的均值递推，与 spp 调度无关 -> 无闪烁 */
        a->mean[i] = (float)(((double)a->mean[i] * N_old + (double)f->mean[i] * f->spp) * inv);
        /* 池化单样本方差；独立样本组均值的方差 = sum(spp*var) / N^2 */
        a->vpool[i] += (float)((double)f->var[i] * f->spp);
        a->mvar[i] = (float)((double)a->vpool[i] / ((double)N_new * N_new));
    }
    a->samples = N_new;
    a->frames++;
    return 0;
}

void rt_frame_passthrough(float* out, const RtFrame* f) {
    memcpy(out, f->mean, (size_t)f->width * f->height * sizeof(float));
}

/* ---------------- 旧接口（仅用于复现闪烁，禁止新代码使用） ---------------- */

/* 采样器侧私有声明：认为自己发送的是“样本和 + 样本计数” */
typedef struct {
    int w, h;
    int frame_no;
    int spp_private;
    const float* sum;     /* 实际放的是均值——语义错位点 1 */
    const float* varmean; /* 实际放的是单样本方差，旧降噪器当均值方差用——错位点 2 */
} LegacySamplerPacket;

/* 降噪器侧私有声明：缓存帧计数从 1 开始，与采样器 0-based 差 1——错位点 3 */
typedef struct {
    float* out;
    float* blend;
    int    frame_count;
} LegacyDenoiser;

static void legacy_init(LegacyDenoiser* d) {
    d->out = calloc(PX, sizeof(float));
    d->blend = calloc(PX, sizeof(float));
    d->frame_count = 1; /* 故意从 1 起：与采样器 frame=0.. 错位 */
}

/* 固定权重时间混合 + 重复归一化：spp 4/8 交替时输出随之明暗交替 */
static void legacy_denoise(LegacyDenoiser* d, const LegacySamplerPacket* p) {
    int spp = spp_of(d->frame_count);           /* 用自己的帧计数推 spp，奇偶可能取反 */
    double alpha = 0.25;
    for (size_t i = 0; i < PX; i++) {
        double incoming = (double)p->sum[i] / (double)spp; /* 均值又除一次 spp */
        if (d->frame_count == 1) d->blend[i] = (float)incoming;
        else d->blend[i] = (float)((1.0 - alpha) * d->blend[i] + alpha * incoming);
        d->out[i] = d->blend[i];
    }
    d->frame_count++;
}

/* ---------------- 指标与产物 ---------------- */

static double image_mean(const float* img) {
    double s = 0;
    for (size_t i = 0; i < PX; i++) s += img[i];
    return s / PX;
}

static double image_mse_to_truth(const float* img) {
    double s = 0;
    for (int y = 0; y < H; y++)
        for (int x = 0; x < W; x++) {
            double e = (double)img[y * W + x] - ground_truth(x, y);
            s += e * e;
        }
    return s / PX;
}

/* 闪烁帧：暖机后帧间整图均值的相对变化超过 FLICKER_EPS */
static int count_flicker(const double* means, int frames, int warmup, int* worst_idx, double* worst_rel) {
    int n = 0;
    *worst_idx = -1; *worst_rel = 0;
    for (int t = warmup + 1; t < frames; t++) {
        double base = fabs(means[t - 1]) > 1e-9 ? fabs(means[t - 1]) : 1e-9;
        double rel = fabs(means[t] - means[t - 1]) / base;
        if (rel > *worst_rel) { *worst_rel = rel; *worst_idx = t; }
        if (rel > FLICKER_EPS) n++;
    }
    return n;
}

static void write_pgm(const char* path, const float* img) {
    FILE* fp = fopen(path, "wb");
    if (!fp) { perror(path); exit(2); }
    fprintf(fp, "P5\n%d %d\n255\n", W, H);
    for (size_t i = 0; i < PX; i++) {
        double v = img[i]; if (v < 0) v = 0; if (v > 1) v = 1;
        unsigned char c = (unsigned char)lround(v * 255.0);
        fwrite(&c, 1, 1, fp);
    }
    fclose(fp);
}

static void write_f32(const char* path, const float* buf, size_t n) {
    FILE* fp = fopen(path, "wb");
    if (!fp) { perror(path); exit(2); }
    unsigned int magic = GOLDEN_MAGIC;
    int w = W, h = H;
    fwrite(&magic, sizeof(magic), 1, fp);
    fwrite(&w, sizeof(int), 1, fp);
    fwrite(&h, sizeof(int), 1, fp);
    fwrite(buf, sizeof(float), n, fp);
    fclose(fp);
}

static int read_f32(const char* path, float* buf, size_t n) {
    FILE* fp = fopen(path, "rb");
    if (!fp) return -1;
    unsigned int magic; int w, h;
    if (fread(&magic, sizeof(magic), 1, fp) != 1 || magic != GOLDEN_MAGIC) { fclose(fp); return -1; }
    if (fread(&w, sizeof(int), 1, fp) != 1 || fread(&h, sizeof(int), 1, fp) != 1 ||
        w != W || h != H) { fclose(fp); return -1; }
    if (fread(buf, sizeof(float), n, fp) != n) { fclose(fp); return -1; }
    fclose(fp);
    return 0;
}

static const char* out_dir(void) {
    const char* d = getenv("DENOISE_OUT");
    return d ? d : "artifacts";
}

static void join_path(char* dst, size_t cap, const char* dir, const char* name) {
    snprintf(dst, cap, "%s/%s", dir, name);
}

/* ---------------- 自检 ---------------- */

static int failures = 0;
static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

int main(int argc, char** argv) {
    int write_golden = (argc > 1 && strcmp(argv[1], "--write-golden") == 0);
    const char* dir = out_dir();
    mkdir(dir, 0755);  /* 已存在则忽略 */

    float* mean = malloc(PX * sizeof(float));
    float* var = malloc(PX * sizeof(float));
    float* sum = malloc(PX * sizeof(float));
    float* pass = malloc(PX * sizeof(float));
    float* noisy0 = malloc(PX * sizeof(float)); /* 冻结第 0 帧噪点，供摘要使用 */
    double legacy_mean[FRAMES], unified_mean[FRAMES];
    double legacy_mse[FRAMES], unified_mse[FRAMES];
    double legacy_pvar[FRAMES], unified_pvar[FRAMES];

    LegacyDenoiser ld; legacy_init(&ld);
    RtAccum acc;
    if (rt_accum_init(&acc, W, H) != 0) { fprintf(stderr, "oom\n"); return 2; }

    rng_state = SEED; /* 两路径共享同一确定性样本序列 */
    for (int t = 0; t < FRAMES; t++) {
        int spp = spp_of(t);
        for (int i = 0; i < PX; i++) {
            double s = 0, s2 = 0;
            int x = i % W, y = i / W;
            for (int k = 0; k < spp; k++) {
                float v = take_sample(x, y);
                s += v; s2 += (double)v * v;
            }
            double m = s / spp;
            double vv = (spp > 1) ? (s2 / spp - m * m) * spp / (spp - 1) : 0.0;
            mean[i] = (float)m;
            var[i] = (float)(vv > 0 ? vv : 0.0);
            sum[i] = (float)m; /* 旧采样器名义上发“和”，实际填均值 */
        }
        if (t == 0) memcpy(noisy0, mean, PX * sizeof(float));

        /* 旧路径 */
        LegacySamplerPacket pkt = { W, H, t, spp, sum, var };
        legacy_denoise(&ld, &pkt);
        legacy_mean[t] = image_mean(ld.out);
        legacy_mse[t] = image_mse_to_truth(ld.out);

        /* 统一路径：RtFrame 契约，字段语义只有 rt_frame.h 一处定义 */
        RtFrame f = { W, H, t, spp, mean, var };
        if (rt_accum_push(&acc, &f) != 0) { fprintf(stderr, "push failed\n"); return 2; }
        unified_mean[t] = image_mean(acc.mean);
        unified_mse[t] = image_mse_to_truth(acc.mean);

        double pv = 0;
        for (int i = 0; i < PX; i++) pv += acc.mvar[i];
        unified_pvar[t] = pv / PX;
        double lp = 0;
        for (int i = 0; i < PX; i++) lp += (double)var[i];
        legacy_pvar[t] = lp / PX; /* 旧路径不做均值方差推导，只透传单样本方差 */
    }

    int lworst, uworst; double lrel, urel;
    int legacy_flicker = count_flicker(legacy_mean, FRAMES, WARMUP, &lworst, &lrel);
    int unified_flicker = count_flicker(unified_mean, FRAMES, WARMUP, &uworst, &urel);

    /* 旁路：关闭统一接口，输出必须与输入逐像素字节一致 */
    RtFrame last = { W, H, FRAMES - 1, spp_of(FRAMES - 1), mean, var };
    rt_frame_passthrough(pass, &last);
    int bypass_identical = memcmp(pass, mean, PX * sizeof(float)) == 0;

    /* 旧路径 2:1 明暗交替的直接证据（暖机后偶帧亮、奇帧暗）。
       EMA 会阻尼振幅，故阈值取 1.05；统一路径同口径比值应接近 1。 */
    double even = 0, odd = 0; int ne = 0, no = 0;
    double ueven = 0, uodd = 0;
    for (int t = WARMUP; t < FRAMES; t++) {
        if (t & 1) { odd += legacy_mean[t]; uodd += unified_mean[t]; no++; }
        else { even += legacy_mean[t]; ueven += unified_mean[t]; ne++; }
    }
    double ratio = (even / ne) / (odd / no);
    double uratio = (ueven / ne) / (uodd / no);

    /* golden：旧路径输出冻结，保证“关闭统一”可逐字节复现 */
    char path[512];
    join_path(path, sizeof(path), dir, "golden_legacy_off.f32");
    int golden_ok;
    if (write_golden) {
        write_f32(path, ld.out, PX);
        golden_ok = 1;
        printf("golden_written=%s\n", path);
    } else {
        float* g = malloc(PX * sizeof(float));
        golden_ok = (read_f32(path, g, PX) == 0 && memcmp(g, ld.out, PX * sizeof(float)) == 0);
        free(g);
    }

    /* 产物：噪点图像（第 0 帧原始输出）、旧路径输出、统一路径输出 */
    join_path(path, sizeof(path), dir, "noisy_frame00.pgm");
    write_pgm(path, noisy0);
    join_path(path, sizeof(path), dir, "legacy_output_frame39.pgm");
    write_pgm(path, ld.out);
    join_path(path, sizeof(path), dir, "unified_output_frame39.pgm");
    write_pgm(path, acc.mean);

    /* 旁路产物：关闭统一接口时的输入与输出，供外部逐字节核验 */
    join_path(path, sizeof(path), dir, "bypass_input.f32");
    write_f32(path, mean, PX);
    join_path(path, sizeof(path), dir, "bypass_output.f32");
    write_f32(path, pass, PX);

    /* 方差统计 CSV：每帧整图均值、对真值 MSE、平均像素方差 */
    join_path(path, sizeof(path), dir, "variance_stats.csv");
    FILE* csv = fopen(path, "w");
    if (!csv) { perror(path); return 2; }
    fprintf(csv, "frame,spp,legacy_mean,unified_mean,legacy_mse,unified_mse,legacy_pixel_var,unified_mean_var\n");
    for (int t = 0; t < FRAMES; t++)
        fprintf(csv, "%d,%d,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f\n",
                t, spp_of(t), legacy_mean[t], unified_mean[t],
                legacy_mse[t], unified_mse[t], legacy_pvar[t], unified_pvar[t]);
    fclose(csv);

    /* ---- check 项 ---- */
    check("contract_size_sane", (int)sizeof(RtFrame) > 24 && (int)sizeof(RtAccum) >= 32);
    check("legacy_flicker_present", legacy_flicker > 0);
    check("unified_flicker_zero", unified_flicker == 0);
    check("legacy_brightness_oscillates", ratio > 1.05);
    check("unified_even_odd_balanced", fabs(uratio - 1.0) < 0.002);
    check("unified_steady_rel_max", urel < FLICKER_EPS);
    check("bypass_pixel_identical", bypass_identical);
    check("unified_variance_converges",
          unified_pvar[FRAMES - 1] < 0.25 * unified_pvar[0]);
    check("unified_mse_below_legacy",
          unified_mse[FRAMES - 1] < legacy_mse[FRAMES - 1]);
    check("legacy_golden_byte_stable", golden_ok);

    printf("summary legacy_flicker_frames=%d unified_flicker_frames=%d "
           "legacy_worst_rel=%.6f unified_worst_rel=%.6f "
           "legacy_even_odd_ratio=%.4f unified_even_odd_ratio=%.4f\n",
           legacy_flicker, unified_flicker, lrel, urel, ratio, uratio);
    printf("summary final_mse legacy=%.6f unified=%.6f final_mean_var unified=%.6f\n",
           legacy_mse[FRAMES - 1], unified_mse[FRAMES - 1], unified_pvar[FRAMES - 1]);
    printf("all_ok=%d\n", failures == 0);

    rt_accum_free(&acc);
    free(ld.out); free(ld.blend);
    free(mean); free(var); free(sum); free(pass); free(noisy0);
    return failures == 0 ? 0 : 1;
}
