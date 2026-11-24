/*
 * CDC 6000 single-precision soft-float (60-bit one's complement).
 *
 * Layout of the low 60 bits of a uint64_t (high 4 bits always clear):
 *
 *   bit 59      sign
 *   bits 58–48  11-bit characteristic (biased exponent)
 *   bits 47–0   48-bit integer coefficient
 *
 * A positive finite word with characteristic X and coefficient C has
 *
 *   value = C * 2^(X - CDC_FP_BIAS - 48)
 *
 * with the binary point to the right of C. CDC_FP_BIAS is 1024
 * (2000 octal), the excess used by the CDC 6000 characteristic.
 * A negative word is the bitwise one's complement of that positive
 * encoding (complement of the whole 60-bit word, not of C alone).
 * +0 is all zeros; -0 is 60 one-bits. This file flushes underflow to
 * a signed zero and maps overflow / invalid to a non-numeric pattern
 * (characteristic 2047, coefficient 0) — not a full hardware
 * indefinite/infinity set.
 *
 * Arithmetic is normalized and truncating (chop toward zero on
 * alignment and on the extra product/quotient bits). No rounding
 * bit, no CDC double-precision word, no COMPASS microcode.
 *
 * Sources for the layout and the excess-1024 / 2000-octal bias:
 *   - Control Data, CDC 6400/6500/6600 Computer Systems Reference
 *     Manual, pub. 60100000 (60-bit FP operand, characteristic bias
 *     2000 octal).
 *   - C. B. Silio, Jr., University of Maryland ENEE notes on CDC 6000
 *     floating point (excess-1024 characteristic, integer coefficient).
 *   - NSSDC CDC 60-bit floating-point data-format notes (1+11+48,
 *     excess 1024), and the same split recorded at
 *     https://60bits.net/ (CDC 6000 floating-point pages).
 *
 * Host double is only the pack/unpack boundary. Every + − × ÷ below
 * runs on the 60-bit word.
 */
#ifndef CDC_FP_H
#define CDC_FP_H

#include <stdint.h>

#define CDC_FP_BITS 60
#define CDC_FP_COEFF_BITS 48
#define CDC_FP_EXP_BITS 11
#define CDC_FP_BIAS 1024

#define CDC_FP_MASK60 ((UINT64_C(1) << CDC_FP_BITS) - UINT64_C(1))
#define CDC_FP_SIGN ((UINT64_C(1) << 59))
#define CDC_FP_COEFF_MASK ((UINT64_C(1) << CDC_FP_COEFF_BITS) - UINT64_C(1))
#define CDC_FP_NORM_BIT (UINT64_C(1) << (CDC_FP_COEFF_BITS - 1))

uint64_t cdc_from_double(double x);
double cdc_to_double(uint64_t w);

/* One's-complement negate of the whole 60-bit word. */
uint64_t cdc_neg(uint64_t w);

uint64_t cdc_add(uint64_t a, uint64_t b);
uint64_t cdc_sub(uint64_t a, uint64_t b);
uint64_t cdc_mul(uint64_t a, uint64_t b);
uint64_t cdc_div(uint64_t a, uint64_t b);

/* -1, 0, +1. +0 and -0 compare equal. */
int cdc_cmp(uint64_t a, uint64_t b);

int cdc_is_zero(uint64_t w);
int cdc_is_indef(uint64_t w);

#endif
