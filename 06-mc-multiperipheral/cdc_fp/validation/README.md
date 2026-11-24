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

`cyber_results.txt` is the FTN5 listing of the same cases, reduced to

```text
CDC_FP_CYBER 1
# id result_oct20
one 17204000000000000000
```

The compare does not ask the host word and the Cyber word to be the
same 60 bits. The host packs `value = coeff * 2^(X - 1024 - 48)`, so
normalized +1 is octal `20014000000000000000` (characteristic 2001).
A 6600 stores the integer-coefficient exponent `E = X - 1072` as a
one's-complement characteristic: `E + 1024` when `E >= 0`, and
`E + 1023` when `E < 0`. For +1, `E = -47` and the characteristic is
1720 octal (`17204000000000000000`), which is the normalized form in
the LLNL LTSS hardware notes. Every finite case in this vector has
`E < 0`, so the stored characteristic is exactly 49 below the host
word and the 48-bit coefficient is the same. `+0` and `-0` match as
bits, including results that chopped to zero. A coefficient mismatch
or any other characteristic fails the compare. There is no ulp slack.

```bash
make -C 06-mc-multiperipheral/cdc_fp check-cyber
```

With the checked-in extract that prints

```text
PASS: cdc_fp_cyber (coeff match, characteristic one's-complement)
```

If `cyber_results.txt` is absent the command prints

```text
PASS: cdc_fp_cyber (skipped, no cyber_results.txt)
```

and exits 0. A present file that disagrees exits 1. The host library
is left on the excess-1024 packing; `cdc_to_double` still round-trips
because pack and unpack share that formula. This file records the
6600 image. It is not a cycle-accurate adder claim.

## FTN5 source

`scripts/emit_ftn.py` reads `cases.txt` and writes
`ftn/cdc_fp_cases.f`. The program is fixed-form Fortran, 72 columns.
It stores each result through `EQUIVALENCE` of a `REAL` and an
`INTEGER`, then prints 20 octal digits with `SHIFT` and `.AND. O"7"`.
Signed zero is `.NOT. 0`. The `(1+ε)-1` powers of one half are `DO`
loops of multiplies, so compile with `OPT=0` and do not let FTN5 fold
them. `make check-ftn` regenerates the file and `cmp`s it, and runs
the listing extractor's self-test (`PASS: cdc_fp_cyber_extract`).
The extractor lowercases ids. CDC display code prints them in
uppercase.

## Running it on DtCyber

This tree does not vendor NOS or the emulator. The checked-in extract
comes from one run of Kevin Jordan's DtCyber
(<https://github.com/kej715/DtCyber>) on the NOS 2.8.7 ready-to-run
image `nos287-full-865` (hobbyist image; not redistributed), deadstarted
as a CYBER 865. The dayfile reported NOS 2.8.7 871/871, FTN5
compilation 0.129 CP seconds, and execution 0.018 CP seconds, then
`STOP`. The deck's `*** CDCFPC COMPLETE` line was not printed; the
58 `CDCVAL` lines were taken from the printer file `LP5xx_C12_E5`.

1. Build DtCyber: `make -f Makefile.linux64 all` from the DtCyber root.
   The emulator opens an X11 display. A headless host can use
   `Xvfb :99 -screen 0 1280x1024x24 -nolisten tcp` and `DISPLAY=:99`.
   The stock `cyber.ini` font is `lucidatypewriter`. If that font is
   missing, XLoadFont aborts. A local `[localfont]` with
   `fontName=Liberation Mono` and `fontType=TrueType` starts. A failure
   to listen on port 175 is a log warning.
2. From `NOS2.8.7`, install the ready-to-run image (`node install rtr`)
   or unzip `nos287-full-865.zip` into that directory so `disks/` and
   `tapes/ds.tap` exist. In `cyber.ini`, point the telnet listener at
   an unprivileged port (`terminals=2323,0x04,16,telnet`) and set
   `cdcnetPrivilegedTcpPortOffset=6600` so the process does not need
   to bind ports below 1024. Keep stdin open: DtCyber treats EOF on
   stdin as shutdown.
3. Submit the batch deck. From any cwd:

```bash
DTCYBER_ROOT=/path/to/DtCyber \
DTCYBER_NOS=/path/to/DtCyber/NOS2.8.7 \
node 06-mc-multiperipheral/cdc_fp/validation/scripts/run_nos_job.js \
  06-mc-multiperipheral/cdc_fp/validation/ftn/cdc_fp_cases.f \
  /tmp/cdcfpc-listing.txt
```

The script deadstarts the `[cyber]` profile, loads the deck on card
reader channel 12 equipment 4 (`CR3447` in `cyber.ini`), and waits for
the dayfile marker. The deck is `USER,GUEST,GUEST` and
`FTN5,L=0,OPT=0,GO`. If that marker does not appear within three
minutes, the script writes `LP5xx_C12_E5` when that file contains
`CDCVAL`. The printer file grows across runs, so extract one job's
block. Printer lines that contain `CDCVAL` are the result words.

4. Turn the listing into the compare file:

```bash
python3 06-mc-multiperipheral/cdc_fp/validation/scripts/extract_cyber_results.py \
  /tmp/cdcfpc-listing.txt \
  -o 06-mc-multiperipheral/cdc_fp/validation/cyber_results.txt
make -C 06-mc-multiperipheral/cdc_fp check-cyber
```

This is a cross-check of the case vector, not a NOS port of VEGAS
or `eemm`.
