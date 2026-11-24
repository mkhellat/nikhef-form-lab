#!/usr/bin/env python3
"""Compare CDC-60 oracle dumps.

Host smoke compares two dumps (typically a fresh oracle run against
validation/host_expected.txt) on the 20-digit octal result word.
--check-known recomputes the excess-1024 packing of the exact cases
and does not trust the C library. --cyber compares a host dump to a
DtCyber FTN5 listing extract; a missing Cyber file is a skip, not a
failure.

Exit 0 prints a PASS line. Exit 1 prints FAIL lines.
"""

import argparse
import sys

MASK60 = (1 << 60) - 1
BIAS = 1024
COEFF_BITS = 48

# glibc-style LCG. Operand = (state mod 48) + 1, four pairs.
LCG_SEED = 20261124
LCG_A = 1103515245
LCG_C = 12345


def oct20(word):
    return f"{word & MASK60:020o}"


def pack_norm(coeff, mag_exp):
    """Pack a positive normalized coefficient. coeff is in [2^47, 2^48)."""
    if coeff < (1 << 47) or coeff >= (1 << COEFF_BITS):
        raise ValueError(f"coeff out of range: {coeff}")
    x = mag_exp + BIAS + COEFF_BITS
    if x <= 0 or x >= 2047:
        raise ValueError(f"characteristic out of range: {x}")
    return ((x << COEFF_BITS) | coeff) & MASK60


def pack_pow2(exp2):
    """Exact 2^exp2 as a CDC single word (exp2 <= 0 for our cases)."""
    return pack_norm(1 << 47, exp2 - 47)


def pack_int(n):
    """Exact positive integer with bit length <= 48."""
    if n == 0:
        return 0
    if n < 0 or n.bit_length() > COEFF_BITS:
        raise ValueError(f"not an exact CDC integer: {n}")
    shift = COEFF_BITS - n.bit_length()
    return pack_norm(n << shift, -shift)


def lcg_pairs(n=4):
    state = LCG_SEED
    pairs = []
    for _ in range(n):
        state = (LCG_A * state + LCG_C) & 0x7FFFFFFF
        a = (state % 48) + 1
        state = (LCG_A * state + LCG_C) & 0x7FFFFFFF
        b = (state % 48) + 1
        pairs.append((a, b))
    return pairs


def parse_dump(path):
    """Return (kind, version, ordered list of dicts)."""
    rows = []
    kind = None
    version = None
    try:
        fh = open(path, "r", encoding="utf-8")
    except FileNotFoundError:
        raise SystemExit(f"FAIL: {path}: no such file")
    with fh:
        for lineno, raw in enumerate(fh, 1):
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if parts[0] in ("CDC_FP_ORACLE", "CDC_FP_CYBER"):
                if kind is not None:
                    raise SystemExit(f"FAIL: {path}:{lineno}: second header")
                if len(parts) != 2 or parts[1] != "1":
                    raise SystemExit(
                        f"FAIL: {path}:{lineno}: want header version 1"
                    )
                kind = parts[0]
                version = parts[1]
                continue
            if len(parts) >= 3 and _is_oct(parts[2]):
                ident, op, word = parts[0], parts[1], parts[2]
            elif len(parts) >= 2 and _is_oct(parts[1]):
                ident, op, word = parts[0], None, parts[1]
            else:
                raise SystemExit(f"FAIL: {path}:{lineno}: bad case line")
            word = word.zfill(20)
            if len(word) != 20:
                raise SystemExit(
                    f"FAIL: {path}:{lineno}: octal word is {len(word)} digits"
                )
            rows.append({"id": ident, "op": op, "oct": word, "line": lineno})
    if kind is None:
        raise SystemExit(f"FAIL: {path}: missing CDC_FP_ORACLE/CYBER 1 header")
    return kind, version, rows


def _is_oct(token):
    return token != "" and all(c in "01234567" for c in token)


def index_rows(rows, path):
    by_id = {}
    for row in rows:
        if row["id"] in by_id:
            raise SystemExit(f"FAIL: {path}: duplicate id {row['id']}")
        by_id[row["id"]] = row
    return by_id


def compare_rows(left_path, left_rows, right_path, right_rows):
    left = index_rows(left_rows, left_path)
    right = index_rows(right_rows, right_path)
    fails = 0
    if set(left) != set(right):
        missing = sorted(set(left) - set(right))
        extra = sorted(set(right) - set(left))
        if missing:
            print(f"FAIL: ids missing from {right_path}: {' '.join(missing)}")
            fails += 1
        if extra:
            print(f"FAIL: ids only in {right_path}: {' '.join(extra)}")
            fails += 1
    for ident in left:
        if ident not in right:
            continue
        if left[ident]["oct"] != right[ident]["oct"]:
            print(
                f"FAIL: {ident} {left_path}={left[ident]['oct']} "
                f"{right_path}={right[ident]['oct']}"
            )
            fails += 1
    return fails


def known_words():
    """Exact result words from the format, not from cdc_fp.c."""
    one = pack_pow2(0)
    two = pack_pow2(1)
    half = pack_pow2(-1)
    quarter = pack_pow2(-2)
    six = pack_int(6)
    m2 = (~pack_int(2)) & MASK60
    e40 = pack_pow2(-40)
    e46 = pack_pow2(-46)
    zplus = 0
    expect = {
        "zplus": zplus,
        "zminus": MASK60,
        "one": one,
        "mone": (~one) & MASK60,
        "two": two,
        "half": half,
        "add11": two,
        "add22": pack_pow2(2),
        "add1m1": zplus,
        "adm2m3": (~pack_int(5)) & MASK60,
        "sub53": two,
        "sub11": zplus,
        "mul23": six,
        "mulm4h": m2,
        "mul05": zplus,
        "div63": two,
        "div12": half,
        "div14": quarter,
        "div42": two,
        "e40": e40,
        # 2^-40 is 128 ulps of 1 (ulp = 2^-47). Aligning it onto
        # mag_exp=-47 shifts by 40, so the addend is 2^7 and the
        # sum stays normalized at mag_exp=-47.
        "sum40": pack_norm((1 << 47) + (1 << 7), -47),
        "can40": e40,
        "e50": pack_pow2(-50),
        "sum50": one,  # 2^-50 is below one ulp of 1
        "can50": zplus,
        "e46": e46,
        # de = 46, addend = 2^47 >> 46 = 2, mag_exp stays -47.
        "sum46": pack_norm((1 << 47) + 2, -47),
        "can46": e46,
        "p7": pack_int(1000003),
        "sq": pack_int(1000003 * 1000003),
    }
    return {ident: oct20(word) for ident, word in expect.items()}


def check_known(path):
    _kind, _version, rows = parse_dump(path)
    by_id = index_rows(rows, path)
    expect = known_words()
    fails = 0
    for ident, word in sorted(expect.items()):
        if ident not in by_id:
            print(f"FAIL: known id {ident} missing from {path}")
            fails += 1
            continue
        got = by_id[ident]["oct"]
        if got != word:
            print(f"FAIL: known {ident} got={got} want={word}")
            fails += 1
    return fails


def check_cases(path):
    """The seeded integers in cases.txt must be the LCG pairs."""
    pairs = lcg_pairs(4)
    want = {}
    for i, (a, b) in enumerate(pairs, 1):
        want[f"r{i}a"] = f"+{a}"
        want[f"r{i}b"] = f"+{b}"
    got = {}
    with open(path, "r", encoding="utf-8") as fh:
        for raw in fh:
            line = raw.split("#", 1)[0].strip()
            if not line:
                continue
            parts = line.split()
            if len(parts) >= 3 and parts[1] == "lit":
                got[parts[0]] = parts[2]
    fails = 0
    for ident, token in sorted(want.items()):
        if got.get(ident) != token:
            print(
                f"FAIL: {path} {ident} lit {got.get(ident)!r} "
                f"want {token} (seed {LCG_SEED})"
            )
            fails += 1
    return fails


def cyber_compare(oracle_path, cyber_path):
    try:
        open(cyber_path, "r", encoding="utf-8").close()
    except FileNotFoundError:
        print(
            "PASS: cdc_fp_cyber "
            "(skipped, no cyber_results.txt)"
        )
        return 0
    _ok, _ov, oracle_rows = parse_dump(oracle_path)
    _ck, _cv, cyber_rows = parse_dump(cyber_path)
    fails = compare_rows(oracle_path, oracle_rows, cyber_path, cyber_rows)
    if fails:
        print(f"FAIL: cdc_fp_cyber ({fails})")
        return 1
    print("PASS: cdc_fp_cyber")
    return 0


def main(argv):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "--check-known",
        metavar="ORACLE",
        help="check exact words against the excess-1024 layout",
    )
    ap.add_argument(
        "--check-cases",
        metavar="CASES",
        help="check seeded lit pairs against the LCG",
    )
    ap.add_argument(
        "--cyber",
        nargs=2,
        metavar=("ORACLE", "CYBER"),
        help="compare host dump to a Cyber extract, or skip if absent",
    )
    ap.add_argument(
        "dumps",
        nargs="*",
        help="two dump files to compare (host smoke)",
    )
    args = ap.parse_args(argv)
    fails = 0
    did = False
    if args.check_cases:
        did = True
        fails += check_cases(args.check_cases)
        if fails == 0:
            print("PASS: cdc_fp_oracle_cases")
    if args.check_known:
        did = True
        n = check_known(args.check_known)
        fails += n
        if n == 0:
            print("PASS: cdc_fp_oracle_known")
    if args.cyber:
        did = True
        fails += cyber_compare(args.cyber[0], args.cyber[1])
    if len(args.dumps) == 2:
        did = True
        left_path, right_path = args.dumps
        _lk, _lv, left_rows = parse_dump(left_path)
        _rk, _rv, right_rows = parse_dump(right_path)
        n = compare_rows(left_path, left_rows, right_path, right_rows)
        fails += n
        if n == 0:
            print("PASS: cdc_fp_oracle")
    elif args.dumps:
        ap.error("host smoke takes exactly two dump files")
    if not did:
        ap.error("nothing to do")
    return 1 if fails else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
