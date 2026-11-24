# Thin VEGAS.  make check  →  PASS: vegas_toy

CC ?= cc
CFLAGS ?= -std=c11 -Wall -Wextra -O2
KINC2 = ../kinc2

.PHONY: check clean

check: test_vegas_toy
	./test_vegas_toy

test_vegas_toy: test_vegas_toy.c vegas.c vegas.h \
		$(KINC2)/ranf.c $(KINC2)/ipow.c $(KINC2)/iipow.c
	$(CC) $(CFLAGS) -I$(KINC2) -Wno-unused-parameter -o $@ \
		test_vegas_toy.c vegas.c \
		$(KINC2)/ranf.c $(KINC2)/ipow.c $(KINC2)/iipow.c -lm

clean:
	rm -f test_vegas_toy
