CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra
# OpenMP parallelizes the path tracer pixel loop; build still works without it.
CFLAGS += -fopenmp
LDFLAGS += -fopenmp -lm

SRC := src/image.c src/texture.c src/tonemap.c src/scene.c src/tri.c \
       src/camera.c src/radiometry.c
LIB := libparity.a
BIN := selftest pathtracer rasterizer diffimg

all: $(BIN)

$(LIB): $(SRC)
	$(CC) $(CFLAGS) -c $(SRC)
	ar rcs $@ *.o

selftest: src/selftest.c $(LIB)
	$(CC) $(CFLAGS) src/selftest.c $(LIB) -o $@ $(LDFLAGS)

pathtracer: src/pathtracer.c $(LIB)
	$(CC) $(CFLAGS) src/pathtracer.c $(LIB) -o $@ $(LDFLAGS)

rasterizer: src/rasterizer.c $(LIB)
	$(CC) $(CFLAGS) src/rasterizer.c $(LIB) -o $@ $(LDFLAGS)

diffimg: src/diffimg.c $(LIB)
	$(CC) $(CFLAGS) src/diffimg.c $(LIB) -o $@ $(LDFLAGS)

# quick smoke test: 64 spp keeps CI fast; full parity run uses the default 4096
test: $(BIN)
	./selftest
	./rasterizer -size 160 120
	./pathtracer -spp 64 -size 160 120
	./diffimg -hf out/rt_hdr.pfm out/rast_hdr.pfm

# full offline render (thousands of samples) then the pixel diff report
parity: $(BIN)
	./rasterizer
	./pathtracer -spp 4096
	./diffimg -hf out/rt_hdr.pfm out/rast_hdr.pfm | tee out/diff_report.txt
	python3 -I tools/histogram.py out/diff_hist.csv out/diff_hist.svg

clean:
	rm -f $(BIN) $(LIB) *.o out/*
