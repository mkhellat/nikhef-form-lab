/*
 * CDC-60 naive |M|^2 and the dressed L⊗L contraction used by
 * kinc2-driver/digit_loss_eemm_gamgam.c.
 *
 * eval_naive_poly_cdc evaluates the 98 monomials with cdc_fp.
 * dressed_naive_cdc / dressed_naive_ld run the same lepton-tensor
 * Hessian. The long-double entry is the IEEE reference; the CDC
 * entry packs each dot product into a 60-bit word and stays there
 * for every arithmetic step of the polynomial. The Hessian combine
 * and the lepton tensor L stay in host double, matching the
 * reference driver, so a mild-kinematics compare isolates the
 * polynomial word length.
 */
#ifndef EVAL_NAIVE_CDC_H
#define EVAL_NAIVE_CDC_H

#include "eval_naive_poly.h"
#include "eval_stable_35.h"

#include <stdint.h>

double eval_naive_poly_cdc(const Dots *d);
uint64_t eval_naive_poly_cdc_w(const Dots *d);

double dressed_naive_cdc(const V4 *q1, const V4 *q2, const V4 *p1,
			 const V4 *p2, const V4 *p6, const V4 *p7, double m);
double dressed_naive_ld(const V4 *q1, const V4 *q2, const V4 *p1,
			const V4 *p2, const V4 *p6, const V4 *p7, double m);

#endif
