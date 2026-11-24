# CDC-60 soft-float.
#   make check        PASS: cdc_fp, PASS: cdc_fp_oracle,
#                     PASS: cdc_fp_cyber or a skip
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

.PHONY: check check-oracle check-cyber oracle clean

check: test_cdc_fp $(ORACLE)
	./test_cdc_fp
	$(MAKE) check-oracle
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

$(ORACLE): validation/oracle_cdc_fp.c cdc_fp.c cdc_fp.h
	$(CC) $(CFLAGS) -I. -o $@ validation/oracle_cdc_fp.c cdc_fp.c -lm

test_cdc_fp: test_cdc_fp.c cdc_fp.c cdc_fp.h
	$(CC) $(CFLAGS) -o $@ test_cdc_fp.c cdc_fp.c -lm

clean:
	rm -f test_cdc_fp $(ORACLE) $(ORACLE_OUT)
