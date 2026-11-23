
/* rt_frame.h: 采样器 <-> 降噪器统一帧接口（唯一契约来源）。
 *
 * 旧实现两边各自声明结构、各自理解字段语义（样本和 vs 均值、单样本方差 vs
 * 均值方差、帧计数错位），导致降噪输出随 spp 调度忽明忽暗。统一后双方只允许
 * 包含本头文件，任何字段语义变更只能改这里。
 *
 * 契约：
 *  - mean[]  像素样本均值（线性辐射度，已按本帧 spp 归一化），不是样本和
 *  - var[]   Bessel 修正的单样本方差（与 sample 模块 Stat.var 同公式）；
 *            累计均值的方差由降噪器按 vpool / N^2 推导，采样器不做除法
 *  - spp     本帧样本数；N 为跨帧累计样本数，由降噪器维护
 *  - frame   帧序号，从 0 开始
 *  - 像素行优先、单通道 float，stride == width，无对齐填充
 */
#ifndef RT_FRAME_H
#define RT_FRAME_H

typedef struct {
    int width;
    int height;
    int frame;              /* 0-based */
    int spp;                /* 本帧样本数 */
    const float* mean;      /* 长度 width*height */
    const float* var;       /* 长度 width*height */
} RtFrame;

typedef struct {
    float* mean;            /* 累计样本均值，长度 w*h */
    float* vpool;           /* 单样本方差的 spp 加权池化值 sum(spp*var) */
    float* mvar;            /* 累计均值方差估计 = vpool / N^2 */
    int    frames;          /* 已累计帧数 */
    long   samples;         /* 累计样本总数 N */
    int    width;
    int    height;
} RtAccum;

int  rt_accum_init(RtAccum* a, int width, int height);
void rt_accum_free(RtAccum* a);

/* 接入一帧，按总样本数加权累计；形状不匹配返回 -1，成功返回 0 */
int  rt_accum_push(RtAccum* a, const RtFrame* f);

/* 关闭统一接口时的旁路：逐像素原样拷贝，调用方输出与输入字节一致 */
void rt_frame_passthrough(float* out, const RtFrame* f);

#endif
