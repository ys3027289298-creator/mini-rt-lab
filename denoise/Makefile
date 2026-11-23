CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra
# RT_UNIFIED=1（默认）走统一接口；=0 仅构建旧路径复现用二进制，产物文件名区分
RT_UNIFIED ?= 1
ifeq ($(RT_UNIFIED),1)
BIN := denoise
else
BIN := denoise_legacy
CFLAGS += -DRT_UNIFIED=0
endif

all: $(BIN)

$(BIN): src/denoise.c src/rt_frame.h
	$(CC) $(CFLAGS) src/denoise.c -o $@ -lm

# 冻结旧路径输出为 golden，供之后的逐字节稳定性校验
golden: src/denoise.c src/rt_frame.h
	$(CC) $(CFLAGS) src/denoise.c -o $(BIN) -lm
	./$(BIN) --write-golden

# 自包含自检：本地 test_out 内先生成 golden 再做字节稳定性校验
test: $(BIN)
	mkdir -p test_out
	DENOISE_OUT=test_out ./$(BIN) --write-golden
	DENOISE_OUT=test_out ./$(BIN)

clean:
	rm -f denoise denoise_legacy
	rm -rf test_out
