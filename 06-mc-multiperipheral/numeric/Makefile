CC ?= cc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra
TARGETS = digit_loss_naive mc_digit_loss

.PHONY: all run mc regen clean

all: $(TARGETS)

digit_loss_naive: digit_loss_naive.c eval_naive_poly.h
	$(CC) $(CFLAGS) -o $@ digit_loss_naive.c -lm

mc_digit_loss: mc_digit_loss.c eval_naive_poly.h
	$(CC) $(CFLAGS) -o $@ mc_digit_loss.c -lm

run: digit_loss_naive
	./digit_loss_naive

mc: mc_digit_loss
	./mc_digit_loss

# Regenerate the 98-term header from sympy (optional; header is committed).
regen:
	cd ../scripts && python3 emit_naive_poly_c.py

clean:
	rm -f $(TARGETS)
