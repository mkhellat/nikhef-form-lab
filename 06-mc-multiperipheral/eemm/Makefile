# PETRA-style VEGAS eemm driver.
#   make check          → PASS: eemm_vegas   (stable Σ M_ij)
#   ./eemm_vegas ieee
#   ./eemm_vegas cdc

CC ?= cc
CFLAGS ?= -std=c11 -Wall -Wextra -O2 -Wno-unused-parameter
KINC2 = ../kinc2
NUM = ../numeric
CDC = ../cdc_fp
ME = ../me
VEGAS = ../vegas

KINC2_SRC = \
	$(KINC2)/orient.c \
	$(KINC2)/pickin.c \
	$(KINC2)/mapt.c \
	$(KINC2)/mapw.c \
	$(KINC2)/mapla.c \
	$(KINC2)/ipow.c \
	$(KINC2)/iipow.c \
	$(KINC2)/ranf.c \
	$(KINC2)/gamgam.c

.PHONY: check clean

check: eemm_vegas
	./eemm_vegas

eemm_vegas: eemm_vegas.c $(VEGAS)/vegas.c $(CDC)/cdc_fp.c \
		$(ME)/eval_naive_cdc.c $(ME)/dressed_naive_cdc.c $(KINC2_SRC)
	$(CC) $(CFLAGS) -I$(KINC2) -I$(NUM) -I$(CDC) -I$(ME) -o $@ \
		eemm_vegas.c $(VEGAS)/vegas.c $(CDC)/cdc_fp.c \
		$(ME)/eval_naive_cdc.c $(ME)/dressed_naive_cdc.c \
		$(KINC2_SRC) -lm

clean:
	rm -f eemm_vegas
