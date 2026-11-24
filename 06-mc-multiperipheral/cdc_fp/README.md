# cdc_fp — CDC 6000 single-precision soft-float

60-bit one's-complement floating point packed in the low 60 bits of a
`uint64_t`. The split is 1 sign bit, an 11-bit characteristic biased by
1024 (2000 octal), and a 48-bit integer coefficient, matching the CDC
6400/6500/6600 reference manual (pub. 60100000), Silio's ENEE notes,
the NSSDC CDC 60-bit format notes, and the layout summarized at
<https://60bits.net/>. A negative word is the bitwise complement of the
whole 60-bit positive encoding. `+0` is all zeros; `-0` is 60 one-bits.

`cdc_add`, `cdc_sub`, `cdc_mul`, and `cdc_div` are normalized and
truncating (chop toward zero). There is no rounding bit, no CDC
double-precision word, and no COMPASS microcode. Pack/unpack against
host `double` is only the boundary; the arithmetic itself stays on the
60-bit word. Overflow becomes a non-numeric pattern (characteristic
2047, coefficient 0). This is not a DtCyber run and not a claim that
the bit patterns match a 6600 adder cycle for cycle.

```bash
make -C 06-mc-multiperipheral/cdc_fp check
```

prints `PASS: cdc_fp`, then `PASS: cdc_fp_oracle_cases`,
`PASS: cdc_fp_oracle_known`, and `PASS: cdc_fp_oracle`. The same run
also compares a fresh dump to `validation/cyber_results.txt` when that
file exists. Without it the line is
`PASS: cdc_fp_cyber (skipped, no cyber_results.txt)`.

`make oracle` writes the case vector to stdout (`CDC_FP_ORACLE 1`,
then `id op oct20 hex15 decimal`). The octal field is the low 60
bits, including one's-complement negatives and signed zero. See
`validation/README.md`.
