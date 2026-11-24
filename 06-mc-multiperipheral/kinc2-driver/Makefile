# Drivers against vendored kinc2/ (π⁰ 2→3 and eemm 2→4).

CC ?= cc
CFLAGS ?= -Wall -g
KINC2 = ../kinc2
NUMERIC = ../numeric
CDC = ../cdc_fp
ME = ../me

KINC2_COMMON = \
	$(KINC2)/orient.c \
	$(KINC2)/pickin.c \
	$(KINC2)/mapt.c \
	$(KINC2)/mapw.c \
	$(KINC2)/mapla.c \
	$(KINC2)/ipow.c \
	$(KINC2)/iipow.c \
	$(KINC2)/ranf.c

.PHONY: all check check-eemm digit-loss-eemm mc-digit-loss-eemm \
	digit-loss-eemm-cdc clean \
	scan_pi0_scaling check_gram_minors mc_xsec_s_scan

all: check

check: check-eemm digit-loss-eemm mc-digit-loss-eemm digit-loss-eemm-cdc

check-eemm: check_eemm_gamgam
	./check_eemm_gamgam

digit-loss-eemm: digit_loss_eemm_gamgam
	./digit_loss_eemm_gamgam

mc-digit-loss-eemm: mc_digit_loss_eemm_gamgam
	./mc_digit_loss_eemm_gamgam

digit-loss-eemm-cdc: digit_loss_eemm_cdc
	./digit_loss_eemm_cdc

check_eemm_gamgam: check_eemm_gamgam.c $(KINC2_COMMON) $(KINC2)/gamgam.c $(KINC2)/peripp.c
	$(CC) $(CFLAGS) -I$(KINC2) -I$(NUMERIC) -o $@ \
		check_eemm_gamgam.c $(KINC2_COMMON) \
		$(KINC2)/gamgam.c $(KINC2)/peripp.c -lm

digit_loss_eemm_gamgam: digit_loss_eemm_gamgam.c $(KINC2_COMMON) $(KINC2)/gamgam.c
	$(CC) $(CFLAGS) -O2 -I$(KINC2) -I$(NUMERIC) -o $@ \
		digit_loss_eemm_gamgam.c $(KINC2_COMMON) \
		$(KINC2)/gamgam.c -lm

mc_digit_loss_eemm_gamgam: mc_digit_loss_eemm_gamgam.c $(KINC2_COMMON) $(KINC2)/gamgam.c
	$(CC) $(CFLAGS) -O2 -I$(KINC2) -I$(NUMERIC) -o $@ \
		mc_digit_loss_eemm_gamgam.c $(KINC2_COMMON) \
		$(KINC2)/gamgam.c -lm

digit_loss_eemm_cdc: digit_loss_eemm_cdc.c $(KINC2_COMMON) $(KINC2)/gamgam.c \
		$(CDC)/cdc_fp.c $(ME)/eval_naive_cdc.c $(ME)/dressed_naive_cdc.c
	$(CC) $(CFLAGS) -O2 -I$(KINC2) -I$(NUMERIC) -I$(CDC) -I$(ME) -o $@ \
		digit_loss_eemm_cdc.c $(KINC2_COMMON) $(KINC2)/gamgam.c \
		$(CDC)/cdc_fp.c $(ME)/eval_naive_cdc.c $(ME)/dressed_naive_cdc.c -lm

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
	rm -f check_eemm_gamgam digit_loss_eemm_gamgam \
		mc_digit_loss_eemm_gamgam digit_loss_eemm_cdc \
		scan_pi0_scaling check_gram_minors mc_xsec_s_scan *.o
