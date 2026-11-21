CC      ?= gcc
CFLAGS  ?= -std=c11 -O1 -g -Wall -Wextra -Wno-unused-function
SAN     := -fsanitize=address,undefined -fno-sanitize-recover=all
BUILD   := build
BINS    := $(BUILD)/extreme_bvh $(BUILD)/extreme_sample $(BUILD)/extreme_material

all: $(BINS)

$(BUILD):
	mkdir -p $(BUILD)

$(BUILD)/extreme_bvh: extreme_bvh.c harness.h ../bvh/src/bvh.c | $(BUILD)
	$(CC) $(CFLAGS) $(SAN) -o $@ extreme_bvh.c -lm

$(BUILD)/extreme_sample: extreme_sample.c harness.h ../sample/src/sample.c | $(BUILD)
	$(CC) $(CFLAGS) $(SAN) -o $@ extreme_sample.c -lm

$(BUILD)/extreme_material: extreme_material.c harness.h ../material/src/material.c | $(BUILD)
	$(CC) $(CFLAGS) $(SAN) -o $@ extreme_material.c -lm

check: all
	@fail=0; \
	for b in $(BINS); do \
	  echo "== $$b =="; \
	  ./$$b || fail=1; \
	done; \
	echo "matrix_all_ok=$$((1 - fail))"; \
	exit $$fail

# 确认矩阵运行后仓库已跟踪文件零改动（样例源码、Makefile 等）
verify-clean:
	@test -z "$$(git -C .. status --porcelain --untracked-files=no)" \
	  && echo "repo_clean=1" \
	  || { echo "repo_modified=1"; git -C .. status --porcelain --untracked-files=no; exit 1; }

clean:
	rm -rf $(BUILD)

.PHONY: all check verify-clean clean
