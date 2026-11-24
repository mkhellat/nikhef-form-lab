# cdc_fp validation harness

Host dump of a fixed CDC-60 case vector, plus a compare against either
that dump or a DtCyber FTN5 extract. The arithmetic under test is
`cdc_add`, `cdc_sub`, `cdc_mul`, and `cdc_div`. This directory does
not claim a cycle-accurate 6600, and it does not run VEGAS.

## Cases

`cases.txt` is the case list. An id is at most 7 characters.

- `lit` packs a decimal token with `cdc_from_double`. `+0` and `-0`
  are signed zeros. The other literals used here are integers or
  dyadics that fit in the 48-bit coefficient.
- `add`, `sub`, `mul`, `div` take two earlier ids.
- `muln A B N` multiplies A by B, N times. The `(1+ε)-1` cases build
  `2^-40`, `2^-46`, and `2^-50` that way, then add to 1 and subtract
  1. `2^-50` is below one ulp of 1, so the add chops it off.
- Four integer pairs come from the LCG
  `s = (1103515245 * s + 12345) mod 2^31` with seed `20261124`.
  Each operand is `(s mod 48) + 1`. They are written as `lit` lines.
  `1/3` and those quotients are truncating divisions, not exact
  decimals.

## Host dump

```bash
make -C 06-mc-multiperipheral/cdc_fp oracle
make -C 06-mc-multiperipheral/cdc_fp check-oracle
```

`oracle` prints:

```text
CDC_FP_ORACLE 1
# id op result_oct20 result_hex15 result_dec
one lit 20014000000000000000 401800000000000 1.0000000000000000e+00
```

`result_oct20` is 20 octal digits (60 bits). `result_hex15` is 15 hex
digits of the same word. The compare uses the octal word only, so a
decimal-format change does not hide a bit difference. `+0` is 20
zeros. `-0` is 20 sevens. Negatives are the one's complement of the
positive word.

`check-oracle` rebuilds the dump and checks three things:

1. The seeded `lit` pairs match the LCG (`PASS: cdc_fp_oracle_cases`).
2. Exact results (`+1`, `2^-40`, `(1+2^-50)-1`, `1000003^2`, …) match
   an independent excess-1024 packing in `compare_cdc_oracle.py`, not
   a second call into `cdc_fp.c` (`PASS: cdc_fp_oracle_known`).
3. The full dump matches `host_expected.txt` on every id's octal word
   (`PASS: cdc_fp_oracle`).

## Cyber extract

When a DtCyber FTN5 run of the same cases exists, put it at
`cyber_results.txt`:

```text
CDC_FP_CYBER 1
# id result_oct20
one 20014000000000000000
```

Ids and 20-digit octal result words must match the host dump exactly,
including the sign of zero. `div` results are held to the same rule.
There is no ulp slack.

```bash
make -C 06-mc-multiperipheral/cdc_fp check-cyber
```

If `cyber_results.txt` is absent the command prints

```text
PASS: cdc_fp_cyber (skipped, no cyber_results.txt)
```

and exits 0. A present file that disagrees exits 1.
