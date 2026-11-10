# FORM checks for the muon-pair cancellation report (Section 5).

FORM ?= form
WARD = mumu_ward.frm
NAIVE = mumu_naive.frm
SCHOUTEN = mumu_schouten_c1.frm
APPC = mumu_schouten_c2_c14.frm
MIJ = mumu_mij.frm

.PHONY: all check ward naive schouten appC mij clean

all: check

ward:
	$(FORM) $(WARD)

naive:
	$(FORM) $(NAIVE)

# Appendix-C toolkit. Match Print assignment lines (leading spaces).
schouten:
	@out=$$(mktemp); \
	$(FORM) $(SCHOUTEN) > $$out 2>&1; \
	ok=1; \
	grep -qE '^[[:space:]]+NnormCheck = 0;' $$out || ok=0; \
	grep -qE '^[[:space:]]+Biline = 0;' $$out || ok=0; \
	grep -qE '^[[:space:]]+SchCheck = 0;' $$out || ok=0; \
	if [ $$ok -eq 1 ]; then \
	  echo "PASS: mumu_schouten_c1.frm"; \
	  rm -f $$out; \
	else \
	  echo "FAIL: expected NnormCheck=0; Biline=0; SchCheck=0;" >&2; \
	  sed -n '1,120p' $$out >&2; \
	  rm -f $$out; \
	  exit 1; \
	fi

# Appendix-C reductions used for M12/M21: (C.4), (C.8), (C.9), (C.13), (C.14).
appC:
	@out=$$(mktemp); \
	$(FORM) $(APPC) > $$out 2>&1; \
	ok=1; \
	grep -qE '^[[:space:]]+DiffC4 = 0;' $$out || ok=0; \
	grep -qE '^[[:space:]]+DiffC8 = 0;' $$out || ok=0; \
	grep -qE '^[[:space:]]+DiffC9 = 0;' $$out || ok=0; \
	grep -qE '^[[:space:]]+DiffC13 = 0;' $$out || ok=0; \
	grep -qE '^[[:space:]]+DiffC14 = 0;' $$out || ok=0; \
	if [ $$ok -eq 1 ]; then \
	  echo "PASS: mumu_schouten_c2_c14.frm"; \
	  rm -f $$out; \
	else \
	  echo "FAIL: expected DiffC4=DiffC8=DiffC9=DiffC13=DiffC14=0;" >&2; \
	  sed -n '1,160p' $$out >&2; \
	  rm -f $$out; \
	  exit 1; \
	fi

# (3.11) M_ij: Levi-Civita form vs closed poly (eval_mij / report).
mij:
	@out=$$(mktemp); \
	$(FORM) $(MIJ) > $$out 2>&1; \
	ok=1; \
	grep -qE '^[[:space:]]+DiffM11 = 0;' $$out || ok=0; \
	grep -qE '^[[:space:]]+DiffM12 = 0;' $$out || ok=0; \
	grep -qE '^[[:space:]]+DiffM21 = 0;' $$out || ok=0; \
	grep -qE '^[[:space:]]+DiffM22 = 0;' $$out || ok=0; \
	if [ $$ok -eq 1 ]; then \
	  echo "PASS: mumu_mij.frm"; \
	  rm -f $$out; \
	else \
	  echo "FAIL: expected DiffM11=DiffM12=DiffM21=DiffM22=0;" >&2; \
	  sed -n '1,200p' $$out >&2; \
	  rm -f $$out; \
	  exit 1; \
	fi

check: $(WARD) $(SCHOUTEN) $(APPC) $(MIJ)
	@out=$$(mktemp); \
	$(FORM) $(WARD) > $$out 2>&1; \
	if grep -q 'Sg = 0;' $$out; then \
	  echo "PASS: mumu_ward.frm -> Sg = 0"; \
	  rm -f $$out; \
	else \
	  echo "FAIL: expected 'Sg = 0;' in FORM output" >&2; \
	  sed -n '1,80p' $$out >&2; \
	  rm -f $$out; \
	  exit 1; \
	fi
	@$(MAKE) --no-print-directory schouten
	@$(MAKE) --no-print-directory appC
	@$(MAKE) --no-print-directory mij

clean:
	rm -f *.log
