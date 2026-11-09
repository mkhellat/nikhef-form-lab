# FORM checks for the muon-pair cancellation report (Section 5).

FORM ?= form
WARD = mumu_ward.frm
NAIVE = mumu_naive.frm
SCHOUTEN = mumu_schouten_c1.frm

.PHONY: all check ward naive schouten clean

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

check: $(WARD) $(SCHOUTEN)
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

clean:
	rm -f *.log
