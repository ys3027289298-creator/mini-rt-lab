PY ?= python3
export PYTHONPATH := $(CURDIR)

.PHONY: all test analyze figs clean

all: test analyze figs

test:
	$(PY) -P test_mse.py

analyze:
	$(PY) -P mse.py

figs:
	$(PY) -P make_figs.py

clean:
	rm -rf __pycache__ results/*.json
