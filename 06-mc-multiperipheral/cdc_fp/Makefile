# CDC-60 soft-float.  make check  →  PASS: cdc_fp

CC ?= cc
CFLAGS ?= -std=c11 -Wall -Wextra -O2

.PHONY: check clean

check: test_cdc_fp
	./test_cdc_fp

test_cdc_fp: test_cdc_fp.c cdc_fp.c cdc_fp.h
	$(CC) $(CFLAGS) -o $@ test_cdc_fp.c cdc_fp.c -lm

clean:
	rm -f test_cdc_fp
