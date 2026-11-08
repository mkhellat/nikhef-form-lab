# FORM checks for the muon-pair cancellation report (Section 5).

FORM ?= form
WARD = mumu_ward.frm
NAIVE = mumu_naive.frm

.PHONY: all check ward naive clean

all: check

ward:
	$(FORM) $(WARD)

naive:
	$(FORM) $(NAIVE)

# Falsifiable gate: squared Ward identity must print Sg = 0.
check: $(WARD)
	@out=$$(mktemp); \
	$(FORM) $(WARD) > $$out 2>&1; \
	if grep -q 'Sg = 0;' $$out; then \
	  echo "PASS: mumu_ward.frm -> Sg = 0"; \
	  rm -f $$out; \
	else \
	  echo "FAIL: expected 'Sg = 0;' in FORM output" >&2; \
	  tail -40 $$out >&2; \
	  rm -f $$out; \
	  exit 1; \
	fi

clean:
	rm -f *.log
