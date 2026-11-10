# Drivers against vendored kinc2/ (π⁰ 2→3 and eemm 2→4).

CC ?= cc
CFLAGS ?= -Wall -g
KINC2 = ../kinc2
NUMERIC = ../numeric

KINC2_COMMON = \
	$(KINC2)/orient.c \
	$(KINC2)/pickin.c \
	$(KINC2)/mapt.c \
	$(KINC2)/mapw.c \
	$(KINC2)/mapla.c \
	$(KINC2)/ipow.c \
	$(KINC2)/iipow.c \
	$(KINC2)/ranf.c

.PHONY: all check check-eemm clean \
	scan_pi0_scaling check_gram_minors mc_xsec_s_scan

all: check-eemm

check: check-eemm

check-eemm: check_eemm_gamgam
	./check_eemm_gamgam

check_eemm_gamgam: check_eemm_gamgam.c $(KINC2_COMMON) $(KINC2)/gamgam.c $(KINC2)/peripp.c
	$(CC) $(CFLAGS) -I$(KINC2) -I$(NUMERIC) -o $@ \
		check_eemm_gamgam.c $(KINC2_COMMON) \
		$(KINC2)/gamgam.c $(KINC2)/peripp.c -lm

scan_pi0_scaling: scan_pi0_scaling.c $(KINC2_COMMON) $(KINC2)/pi0.c
	$(CC) $(CFLAGS) -I$(KINC2) -o $@ \
		scan_pi0_scaling.c $(KINC2_COMMON) $(KINC2)/pi0.c -lm

check_gram_minors: check_gram_minors.c $(KINC2_COMMON) $(KINC2)/pi0.c
	$(CC) $(CFLAGS) -I$(KINC2) -o $@ \
		check_gram_minors.c $(KINC2_COMMON) $(KINC2)/pi0.c -lm

mc_xsec_s_scan: mc_xsec_s_scan.c $(KINC2_COMMON) $(KINC2)/pi0.c
	$(CC) $(CFLAGS) -I$(KINC2) -o $@ \
		mc_xsec_s_scan.c $(KINC2_COMMON) $(KINC2)/pi0.c -lm

clean:
	rm -f check_eemm_gamgam scan_pi0_scaling check_gram_minors mc_xsec_s_scan *.o
