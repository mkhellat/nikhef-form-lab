# CDC-naive |M|^2 adapter.  make check  →  PASS: compare_naive_cdc

CC ?= cc
CFLAGS ?= -std=c11 -Wall -Wextra -O2
CDC = ../cdc_fp
NUM = ../numeric

.PHONY: check clean

check: compare_naive_cdc
	./compare_naive_cdc

compare_naive_cdc: compare_naive_cdc.c eval_naive_cdc.c dressed_naive_cdc.c \
		$(CDC)/cdc_fp.c eval_naive_cdc.h
	$(CC) $(CFLAGS) -I$(CDC) -I$(NUM) -o $@ \
		compare_naive_cdc.c eval_naive_cdc.c dressed_naive_cdc.c \
		$(CDC)/cdc_fp.c -lm

clean:
	rm -f compare_naive_cdc
