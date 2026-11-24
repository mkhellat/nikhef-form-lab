#include "cdc_fp.h"

#include <math.h>
#include <stdint.h>

/*
 * Internal magnitude: value = ±coeff * 2^mag_exp, coeff an integer.
 * Normalized nonzero coeffs sit in [2^47, 2^48).
 * sign 0 = positive, 1 = negative. indef is the overflow pattern.
 */
typedef struct {
	int sign;
	int indef;
	int mag_exp;
	uint64_t coeff;
} CdcMag;

static uint64_t cdc_indef(int sign)
{
	uint64_t w = (UINT64_C(0x7FF) << CDC_FP_COEFF_BITS);
	if (sign)
		w = (~w) & CDC_FP_MASK60;
	return w;
}

static CdcMag cdc_unpack(uint64_t w)
{
	CdcMag m;
	int X;

	m.sign = 0;
	m.indef = 0;
	m.mag_exp = 0;
	m.coeff = 0;
	w &= CDC_FP_MASK60;
	if (w == 0)
		return m;
	if (w == CDC_FP_MASK60) {
		m.sign = 1;
		return m;
	}
	if (w & CDC_FP_SIGN) {
		m.sign = 1;
		w = (~w) & CDC_FP_MASK60;
	}
	X = (int)((w >> CDC_FP_COEFF_BITS) & 0x7FF);
	m.coeff = w & CDC_FP_COEFF_MASK;
	if (X >= 2047 || (X == 0 && m.coeff == 0)) {
		m.indef = (X >= 2047);
		m.coeff = 0;
		return m;
	}
	m.mag_exp = X - CDC_FP_BIAS - CDC_FP_COEFF_BITS;
	return m;
}

static uint64_t cdc_pack(CdcMag m)
{
	int X;

	if (m.indef)
		return cdc_indef(m.sign);
	if (m.coeff == 0)
		return m.sign ? CDC_FP_MASK60 : UINT64_C(0);

	while (m.coeff >= (UINT64_C(1) << CDC_FP_COEFF_BITS)) {
		m.coeff >>= 1;
		m.mag_exp += 1;
	}
	while (m.coeff > 0 && m.coeff < CDC_FP_NORM_BIT) {
		if (m.mag_exp <= (1 - CDC_FP_BIAS - CDC_FP_COEFF_BITS))
			return m.sign ? CDC_FP_MASK60 : UINT64_C(0);
		m.coeff <<= 1;
		m.mag_exp -= 1;
	}
	X = m.mag_exp + CDC_FP_BIAS + CDC_FP_COEFF_BITS;
	if (X <= 0)
		return m.sign ? CDC_FP_MASK60 : UINT64_C(0);
	if (X >= 2047)
		return cdc_indef(m.sign);
	{
		uint64_t w = ((uint64_t)X << CDC_FP_COEFF_BITS) |
			     (m.coeff & CDC_FP_COEFF_MASK);
		if (m.sign)
			w = (~w) & CDC_FP_MASK60;
		return w;
	}
}

uint64_t cdc_neg(uint64_t w)
{
	return (~w) & CDC_FP_MASK60;
}

int cdc_is_indef(uint64_t w)
{
	w &= CDC_FP_MASK60;
	if (w & CDC_FP_SIGN)
		w = (~w) & CDC_FP_MASK60;
	return ((w >> CDC_FP_COEFF_BITS) & 0x7FF) == 0x7FF;
}

int cdc_is_zero(uint64_t w)
{
	CdcMag m;

	if (cdc_is_indef(w))
		return 0;
	m = cdc_unpack(w);
	return m.coeff == 0 && !m.indef;
}

uint64_t cdc_from_double(double x)
{
	CdcMag m;
	int exp;
	double frac;
	double scaled;

	m.sign = 0;
	m.indef = 0;
	m.mag_exp = 0;
	m.coeff = 0;
	if (isnan(x))
		return cdc_indef(0);
	if (x == 0.0)
		return signbit(x) ? CDC_FP_MASK60 : UINT64_C(0);
	if (x < 0.0) {
		m.sign = 1;
		x = -x;
	}
	if (isinf(x))
		return cdc_indef(m.sign);
	frac = frexp(x, &exp);
	scaled = frac * 0x1p48;
	if (scaled >= 0x1p48) {
		m.coeff = UINT64_C(1) << (CDC_FP_COEFF_BITS - 1);
		m.mag_exp = exp - (CDC_FP_COEFF_BITS - 1);
	} else {
		m.coeff = (uint64_t)scaled;
		m.mag_exp = exp - CDC_FP_COEFF_BITS;
	}
	return cdc_pack(m);
}

double cdc_to_double(uint64_t w)
{
	CdcMag m = cdc_unpack(w);
	double v;

	if (m.indef)
		return m.sign ? -HUGE_VAL : HUGE_VAL;
	if (m.coeff == 0)
		return m.sign ? -0.0 : 0.0;
	v = ldexp((double)m.coeff, m.mag_exp);
	return m.sign ? -v : v;
}

static uint64_t cdc_add_mag(CdcMag a, CdcMag b)
{
	CdcMag r;
	int de;
	uint64_t cb;

	if (a.indef || b.indef) {
		r = a;
		r.indef = 1;
		r.coeff = 0;
		return cdc_pack(r);
	}
	if (a.coeff == 0)
		return cdc_pack(b);
	if (b.coeff == 0)
		return cdc_pack(a);
	if (a.mag_exp < b.mag_exp) {
		CdcMag t = a;
		a = b;
		b = t;
	}
	de = a.mag_exp - b.mag_exp;
	if (de >= CDC_FP_COEFF_BITS)
		cb = 0;
	else
		cb = b.coeff >> de;
	r.sign = a.sign;
	r.indef = 0;
	r.mag_exp = a.mag_exp;
	if (a.sign == b.sign) {
		r.coeff = a.coeff + cb;
	} else if (a.coeff >= cb) {
		r.coeff = a.coeff - cb;
		if (r.coeff == 0)
			r.sign = 0;
	} else {
		r.coeff = cb - a.coeff;
		r.sign = b.sign;
	}
	return cdc_pack(r);
}

uint64_t cdc_add(uint64_t a, uint64_t b)
{
	return cdc_add_mag(cdc_unpack(a), cdc_unpack(b));
}

uint64_t cdc_sub(uint64_t a, uint64_t b)
{
	return cdc_add(a, cdc_neg(b));
}

static void mul48(uint64_t a, uint64_t b, uint64_t *hi, uint64_t *lo)
{
	uint64_t a0 = a & 0xFFFFFFFFULL;
	uint64_t a1 = a >> 32;
	uint64_t b0 = b & 0xFFFFFFFFULL;
	uint64_t b1 = b >> 32;
	uint64_t p0 = a0 * b0;
	uint64_t p1 = a1 * b0;
	uint64_t p2 = a0 * b1;
	uint64_t p3 = a1 * b1;
	uint64_t mid = (p0 >> 32) + (p1 & 0xFFFFFFFFULL) + (p2 & 0xFFFFFFFFULL);

	*lo = (p0 & 0xFFFFFFFFULL) | (mid << 32);
	*hi = p3 + (p1 >> 32) + (p2 >> 32) + (mid >> 32);
}

/* Bit length of a 96-bit product stored in hi:lo (hi is bits 64+). */
static int bitlen96(uint64_t hi, uint64_t lo)
{
	int n;

	if (hi) {
		n = 64;
		while (hi > 1) {
			hi >>= 1;
			n++;
		}
		return n;
	}
	if (!lo)
		return 0;
	n = 0;
	while (lo > 1) {
		lo >>= 1;
		n++;
	}
	return n;
}

static uint64_t shr96(uint64_t hi, uint64_t lo, int s, uint64_t *out_lo)
{
	if (s <= 0) {
		*out_lo = lo;
		return hi;
	}
	if (s >= 64) {
		*out_lo = hi >> (s - 64);
		return 0;
	}
	*out_lo = (lo >> s) | (hi << (64 - s));
	return hi >> s;
}

uint64_t cdc_mul(uint64_t a, uint64_t b)
{
	CdcMag aa = cdc_unpack(a);
	CdcMag bb = cdc_unpack(b);
	CdcMag r;
	uint64_t hi, lo;
	int bl, S;

	r.sign = aa.sign ^ bb.sign;
	r.indef = 0;
	r.mag_exp = 0;
	r.coeff = 0;
	if (aa.indef || bb.indef) {
		r.indef = 1;
		return cdc_pack(r);
	}
	if (aa.coeff == 0 || bb.coeff == 0)
		return cdc_pack(r);
	mul48(aa.coeff, bb.coeff, &hi, &lo);
	bl = bitlen96(hi, lo);
	/* Shift so the leading bit sits at position 48; cdc_pack then
	   chops one more bit into the normalized [2^47, 2^48) window
	   and bumps mag_exp. Net shift is (bitlen - 47). */
	S = bl - CDC_FP_COEFF_BITS;
	if (S < 0)
		S = 0;
	(void)shr96(hi, lo, S, &r.coeff);
	r.mag_exp = aa.mag_exp + bb.mag_exp + S;
	return cdc_pack(r);
}

/*
 * Quotient of a 96-bit dividend (hi:lo) by a 48-bit divisor.
 * Caller keeps hi < div so the quotient fits in 64 bits.
 * Chop (no rounding).
 */
static uint64_t div96by48(uint64_t hi, uint64_t lo, uint64_t div)
{
	uint64_t q = 0;
	uint64_t r = hi;
	int i;

	for (i = 0; i < 64; i++) {
		uint64_t bit = lo >> 63;
		int top;

		lo <<= 1;
		q <<= 1;
		top = (int)(r >> 63);
		r = (r << 1) | bit;
		if (top || r >= div) {
			r -= div;
			q |= 1ULL;
		}
	}
	return q;
}

uint64_t cdc_div(uint64_t a, uint64_t b)
{
	CdcMag aa = cdc_unpack(a);
	CdcMag bb = cdc_unpack(b);
	CdcMag r;
	uint64_t q;
	int s;

	r.sign = aa.sign ^ bb.sign;
	r.indef = 0;
	r.mag_exp = 0;
	r.coeff = 0;
	if (aa.indef || bb.indef || bb.coeff == 0) {
		r.indef = 1;
		return cdc_pack(r);
	}
	if (aa.coeff == 0)
		return cdc_pack(r);
	/* (coeff_a << 48) / coeff_b, both coeffs < 2^48. */
	q = div96by48(aa.coeff >> 16, aa.coeff << 48, bb.coeff);
	s = 0;
	if (q >= (UINT64_C(1) << CDC_FP_COEFF_BITS)) {
		q >>= 1;
		s = 1;
	}
	r.coeff = q;
	r.mag_exp = aa.mag_exp - bb.mag_exp - CDC_FP_COEFF_BITS + s;
	return cdc_pack(r);
}

int cdc_cmp(uint64_t a, uint64_t b)
{
	CdcMag aa, bb;

	if (cdc_is_indef(a) || cdc_is_indef(b))
		return 0;
	aa = cdc_unpack(a);
	bb = cdc_unpack(b);
	if (aa.coeff == 0 && bb.coeff == 0)
		return 0;
	if (aa.coeff == 0)
		return bb.sign ? 1 : -1;
	if (bb.coeff == 0)
		return aa.sign ? -1 : 1;
	if (aa.sign != bb.sign)
		return aa.sign ? -1 : 1;
	if (aa.mag_exp != bb.mag_exp)
		return ((aa.mag_exp > bb.mag_exp) ^ aa.sign) ? 1 : -1;
	if (aa.coeff == bb.coeff)
		return 0;
	return ((aa.coeff > bb.coeff) ^ aa.sign) ? 1 : -1;
}
