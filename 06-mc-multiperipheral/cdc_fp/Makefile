# CDC-60 soft-float.
#   make check        PASS: cdc_fp, PASS: cdc_fp_oracle*,
#                     PASS: cdc_fp_cyber_extract,
#                     PASS: cdc_fp_cyber (one's-complement characteristic,
#                     or a skip when cyber_results.txt is absent)
#   make oracle       dump validation/cases.txt to stdout
#   make check-cyber  compare a fresh dump to validation/cyber_results.txt

CC ?= cc
CFLAGS ?= -std=c11 -Wall -Wextra -O2
PYTHON ?= python3

ORACLE = validation/oracle_cdc_fp
ORACLE_OUT = validation/.oracle_out.txt
CASES = validation/cases.txt
EXPECT = validation/host_expected.txt
CYBER = validation/cyber_results.txt

.PHONY: check check-oracle check-cyber check-ftn oracle clean

FTN = validation/ftn/cdc_fp_cases.f
FTN_GEN = validation/.ftn_gen.f

check: test_cdc_fp $(ORACLE)
	./test_cdc_fp
	$(MAKE) check-oracle
	$(MAKE) check-ftn
	$(MAKE) check-cyber

oracle: $(ORACLE)
	./$(ORACLE) $(CASES)

check-oracle: $(ORACLE)
	./$(ORACLE) $(CASES) > $(ORACLE_OUT)
	$(PYTHON) validation/compare_cdc_oracle.py --check-cases $(CASES)
	$(PYTHON) validation/compare_cdc_oracle.py --check-known $(ORACLE_OUT)
	$(PYTHON) validation/compare_cdc_oracle.py $(ORACLE_OUT) $(EXPECT)

check-cyber: $(ORACLE)
	./$(ORACLE) $(CASES) > $(ORACLE_OUT)
	$(PYTHON) validation/compare_cdc_oracle.py --cyber $(ORACLE_OUT) $(CYBER)

check-ftn:
	$(PYTHON) validation/scripts/emit_ftn.py $(CASES) $(FTN_GEN)
	cmp $(FTN_GEN) $(FTN)
	$(PYTHON) validation/scripts/extract_cyber_results.py --self-test

$(ORACLE): validation/oracle_cdc_fp.c cdc_fp.c cdc_fp.h
	$(CC) $(CFLAGS) -I. -o $@ validation/oracle_cdc_fp.c cdc_fp.c -lm

test_cdc_fp: test_cdc_fp.c cdc_fp.c cdc_fp.h
	$(CC) $(CFLAGS) -o $@ test_cdc_fp.c cdc_fp.c -lm

clean:
	rm -f test_cdc_fp $(ORACLE) $(ORACLE_OUT) $(FTN_GEN)
