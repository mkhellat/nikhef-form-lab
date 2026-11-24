/*
 * Golden patterns and cancellation stress for cdc_fp.
 * Pass line: PASS: cdc_fp
 */
#include "cdc_fp.h"

#include <math.h>
#include <stdio.h>
#include <stdint.h>

static int g_fail;

static void expect(int cond, const char *msg)
{
	if (!cond) {
		fprintf(stderr, "FAIL: %s\n", msg);
		g_fail++;
	}
}

static int near(double a, double b, double rel)
{
	double scale = fabs(a) > fabs(b) ? fabs(a) : fabs(b);
	if (scale < 1e-300)
		return fabs(a - b) < 1e-300;
	return fabs(a - b) / scale <= rel;
}

static void expect_d(uint64_t w, double want, const char *msg)
{
	double got = cdc_to_double(w);
	if (!near(got, want, 0x1p-47)) {
		fprintf(stderr, "FAIL: %s got=%.17g want=%.17g\n", msg, got,
			want);
		g_fail++;
	}
}

int main(void)
{
	uint64_t p1, m1, one, two, three, six;
	uint64_t w;
	double x;

	expect(cdc_from_double(0.0) == 0, "+0 word is 0");
	expect(cdc_from_double(-0.0) == CDC_FP_MASK60, "-0 word is 60 ones");
	expect(cdc_is_zero(0), "+0 is zero");
	expect(cdc_is_zero(CDC_FP_MASK60), "-0 is zero");
	expect(cdc_to_double(0) == 0.0, "unpack +0");
	expect(signbit(cdc_to_double(CDC_FP_MASK60)), "unpack -0 sign");
	expect(cdc_cmp(0, CDC_FP_MASK60) == 0, "+0 == -0");

	one = cdc_from_double(1.0);
	/* +1: characteristic 1025, coefficient 2^47. */
	p1 = (UINT64_C(1025) << 48) | (UINT64_C(1) << 47);
	expect(one == p1, "+1 bit pattern excess-1024");
	expect(cdc_to_double(one) == 1.0, "unpack +1");
	m1 = cdc_neg(one);
	expect(m1 == ((~p1) & CDC_FP_MASK60), "-1 is complement of +1");
	expect(cdc_to_double(m1) == -1.0, "unpack -1");
	expect(cdc_from_double(-1.0) == m1, "pack -1");

	two = cdc_from_double(2.0);
	three = cdc_from_double(3.0);
	expect_d(two, 2.0, "pack 2");
	expect_d(three, 3.0, "pack 3");
	expect_d(cdc_from_double(0.5), 0.5, "pack 1/2");
	expect_d(cdc_from_double(-7.0), -7.0, "pack -7");

	expect_d(cdc_add(one, one), 2.0, "1+1");
	expect_d(cdc_add(two, two), 4.0, "2+2");
	expect_d(cdc_add(one, m1), 0.0, "1+(-1)");
	expect_d(cdc_add(cdc_from_double(-2.0), cdc_from_double(-3.0)), -5.0,
		 "-2+-3");
	expect_d(cdc_sub(cdc_from_double(5.0), three), 2.0, "5-3");
	expect(cdc_is_zero(cdc_sub(one, one)), "1-1 is zero");

	six = cdc_mul(two, three);
	expect_d(six, 6.0, "2*3");
	expect_d(cdc_mul(cdc_from_double(-4.0), cdc_from_double(0.5)), -2.0,
		 "(-4)*0.5");
	expect(cdc_is_zero(cdc_mul(cdc_from_double(0.0), cdc_from_double(5.0))),
	       "0*5");
	expect_d(cdc_mul(three, three), 9.0, "3*3");

	expect_d(cdc_div(six, three), 2.0, "6/3");
	expect_d(cdc_div(one, two), 0.5, "1/2");
	expect_d(cdc_div(one, cdc_from_double(4.0)), 0.25, "1/4");
	expect_d(cdc_div(cdc_from_double(4.0), two), 2.0, "4/2");
	expect(cdc_is_indef(cdc_div(one, cdc_from_double(0.0))), "1/0 indef");

	/* 2^-40 sits well inside the 48-bit coefficient of 1 (ulp = 2^-47). */
	w = cdc_sub(cdc_add(one, cdc_from_double(0x1p-40)), one);
	expect_d(w, 0x1p-40, "(1+2^-40)-1");

	/* 2^-50 is below one ulp of 1, so the add chops it away. */
	w = cdc_add(one, cdc_from_double(0x1p-50));
	expect(w == one, "1+2^-50 chops to 1");
	expect(cdc_is_zero(cdc_sub(w, one)), "(1+2^-50)-1 is 0");

	/* 2^-46 is two ulps of 1 and must survive. */
	w = cdc_sub(cdc_add(one, cdc_from_double(0x1p-46)), one);
	expect_d(w, 0x1p-46, "(1+2^-46)-1");

	expect(cdc_cmp(one, two) < 0, "1<2");
	expect(cdc_cmp(two, one) > 0, "2>1");
	expect(cdc_cmp(cdc_from_double(-3.0), one) < 0, "-3<1");
	expect(cdc_cmp(cdc_from_double(-1.0), cdc_from_double(-2.0)) > 0,
	       "-1>-2");
	expect(cdc_cmp(three, three) == 0, "3==3");

	x = 12345.6789012345;
	expect(near(cdc_to_double(cdc_from_double(x)), x, 0x1p-46),
	       "roundtrip ~14 digits");
	x = -0.875 * 0x1p-10;
	expect(cdc_to_double(cdc_from_double(x)) == x, "dyadic -0.875*2^-10");

	if (g_fail) {
		fprintf(stderr, "FAIL: cdc_fp (%d)\n", g_fail);
		return 1;
	}
	printf("PASS: cdc_fp\n");
	return 0;
}
