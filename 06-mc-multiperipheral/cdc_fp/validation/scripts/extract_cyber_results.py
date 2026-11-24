#!/usr/bin/env python3
"""Pull CDCVAL lines out of an FTN5 terminal or printer log.

Writes a CDC_FP_CYBER 1 dump. 20I1 output may be contiguous or split
by spaces; digits are joined. Ids are the first token after CDCVAL.
CDC display code prints those ids in uppercase, so they are stored
in lowercase to match cases.txt.
"""

import argparse
import sys


def extract(text):
    rows = []
    seen = set()
    for lineno, line in enumerate(text.splitlines(), 1):
        if "CDCVAL" not in line:
            continue
        rest = line.split("CDCVAL", 1)[1].strip()
        if not rest:
            raise SystemExit(f"FAIL: line {lineno}: CDCVAL with no payload")
        parts = rest.split()
        ident = parts[0].lower()
        if not ident.isalnum() or not 1 <= len(ident) <= 7:
            raise SystemExit(f"FAIL: line {lineno}: bad id {parts[0]!r}")
        digits = "".join(ch for ch in "".join(parts[1:]) if ch in "01234567")
        if len(digits) != 20:
            raise SystemExit(
                f"FAIL: line {lineno}: {ident} has {len(digits)} octal digits"
            )
        if ident in seen:
            raise SystemExit(f"FAIL: duplicate CDCVAL id {ident}")
        seen.add(ident)
        rows.append((ident, digits))
    if not rows:
        raise SystemExit("FAIL: no CDCVAL lines")
    return rows


def render(rows):
    out = ["CDC_FP_CYBER 1", "# id result_oct20"]
    for ident, digits in rows:
        out.append(f"{ident} {digits}")
    return "\n".join(out) + "\n"


def self_test():
    sample = "\n".join(
        [
            " junk",
            " CDCVAL zplus   00000000000000000000",
            "CDCVAL ONE 20014000000000000000",
            " CDCVAL half  2 0 0 0 4 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0",
        ]
    )
    rows = extract(sample)
    if rows != [
        ("zplus", "00000000000000000000"),
        ("one", "20014000000000000000"),
        ("half", "20004000000000000000"),
    ]:
        raise SystemExit(f"FAIL: extract self-test {rows}")
    text = render(rows)
    if not text.startswith("CDC_FP_CYBER 1\n"):
        raise SystemExit("FAIL: header")
    print("PASS: cdc_fp_cyber_extract")


def main(argv):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--self-test", action="store_true")
    ap.add_argument("log", nargs="?")
    ap.add_argument("-o", "--output")
    args = ap.parse_args(argv)
    if args.self_test:
        self_test()
        if args.log is None:
            return 0
    if not args.log:
        ap.error("log path required")
    with open(args.log, "r", encoding="utf-8", errors="replace") as fh:
        text = render(extract(fh.read()))
    if args.output:
        with open(args.output, "w", encoding="utf-8", newline="\n") as fh:
            fh.write(text)
    else:
        sys.stdout.write(text)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
