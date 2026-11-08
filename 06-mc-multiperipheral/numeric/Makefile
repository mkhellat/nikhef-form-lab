CC ?= cc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra
TARGET = digit_loss_naive

.PHONY: all run regen clean

all: $(TARGET)

$(TARGET): digit_loss_naive.c eval_naive_poly.h
	$(CC) $(CFLAGS) -o $@ digit_loss_naive.c -lm

run: $(TARGET)
	./$(TARGET)

# Regenerate the 98-term header from sympy (optional; header is committed).
regen:
	cd ../scripts && python3 emit_naive_poly_c.py

clean:
	rm -f $(TARGET)
