CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra

all: pipeline

pipeline: src/pipeline.c
	$(CC) $(CFLAGS) $< -o $@ -lm

# 缺陷复现版：累积/降噪两处归一化各多除一个 2（共 4 倍能量损失）
pipeline-buggy: src/pipeline.c
	$(CC) $(CFLAGS) -DENERGY_BUGS $< -o $@ -lm

# outdir [cut]: 默认跑全部 5 个 cut，输出 PFM/PPM/亮度文本到 outdir
run: pipeline
	mkdir -p $(OUT)
	./pipeline $(OUT) $(CUT)

test: pipeline
	./pipeline /tmp/pipeline_selftest

clean:
	rm -f pipeline pipeline-buggy
