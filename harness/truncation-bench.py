#!/usr/bin/env python3
"""Self-supervised truncation benchmark for ram-extend.

Take a curated seed alignment (Stockholm, member names
<seqid>:<start>-<end>_<orient>, 1-based inclusive), pull every member's
genomic range in by T bp on each side, run ram-extend on the truncated
ranges, and score how much of the removed sequence each member's extension
recovered — and how far it over-extended past the original (curated) ends.

This scores the extension engine independent of seeding: the truth is the
curated alignment itself.

Metrics per member and per family:
  recovery  = recovered bp / T, per side (1.0 = perfect)
  overshoot = bp extended beyond the curated end, per side

Usage:
  truncation-bench.py <in.stk> <genome.2bit> <T> [workdir]
      [-- extra ram-extend args]

Defaults mirror extend-stk.pl: -L 20000 -bandwidth 14 -matrix 20p43g
-minimprovement 27 (override after --).
"""
import os
import re
import subprocess
import sys
import tempfile

RAM_EXTEND = os.environ.get(
    "RAM_EXTEND",
    os.path.join(os.path.dirname(__file__), "..", "target", "release", "ram-extend"),
)
DEFAULT_ARGS = ["-L", "20000", "-bandwidth", "14", "-matrix", "20p43g",
                "-minimprovement", "27"]
# Skip members whose truncated core would be shorter (env-overridable for
# short-element families).
MIN_CORE = int(os.environ.get("TRUNC_MIN_CORE", "50"))


def parse_stk_members(stk_path):
    members = []
    with open(stk_path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line or line.startswith("#") or line.startswith("//"):
                continue
            name = line.split(None, 1)[0]
            m = re.match(r"^(\S+):(\d+)-(\d+)_([+-])$", name)
            if not m:
                sys.exit(f"Unparseable sequence name: {name}")
            members.append((m.group(1), int(m.group(2)), int(m.group(3)), m.group(4)))
    return members


def parse_outtsv(tsv_path):
    """Map (seqid, anchor_lo, anchor_hi) -> (ext_start, ext_end, orient)."""
    hits = {}
    with open(tsv_path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 5:
                continue
            seqid, ext_s, ext_e, orient, details = f[0], int(f[1]), int(f[2]), f[3], f[4]
            m = re.search(r"anchor_range=(\d+)-(\d+)", details)
            if not m:
                continue
            a1, a2 = int(m.group(1)), int(m.group(2))
            hits[(seqid, min(a1, a2), max(a1, a2))] = (ext_s, ext_e, orient)
    return hits


def main():
    args = sys.argv[1:]
    extra = DEFAULT_ARGS
    if "--" in args:
        i = args.index("--")
        args, extra = args[:i], args[i + 1:]
    if len(args) < 3:
        sys.exit(__doc__)
    stk_path, twobit, t = args[0], args[1], int(args[2])
    work = args[3] if len(args) > 3 else tempfile.mkdtemp(prefix="truncbench.")
    os.makedirs(work, exist_ok=True)

    members = parse_stk_members(stk_path)
    ranges_path = os.path.join(work, "trunc.tsv")
    truth = []  # (seqid, true_start, true_end, orient, anchor_key)
    with open(ranges_path, "w") as out:
        for seqid, start, end, orient in members:
            if end - start + 1 < 2 * t + MIN_CORE:
                continue
            ts, te = start + t, end - t  # truncated, 1-based inclusive
            # ranges TSV wants 0-based half-open start
            out.write(f"{seqid}\t{ts - 1}\t{te}\t1\t1\t{orient}\n")
            truth.append((seqid, start, end, orient, (seqid, ts, te)))
    if not truth:
        sys.exit(f"no members long enough for T={t}")

    out_tsv = os.path.join(work, "ext.tsv")
    log = os.path.join(work, "ram-extend.log")
    cmd = [RAM_EXTEND, "-twobit", twobit, "-ranges", ranges_path,
           "-outtsv", out_tsv] + extra
    with open(log, "w") as lg:
        subprocess.run(cmd, stdout=lg, stderr=subprocess.STDOUT, check=True)

    hits = parse_outtsv(out_tsv)
    rows = []
    for seqid, true_s, true_e, orient, key in truth:
        if key not in hits:
            rows.append((seqid, true_s, true_e, orient, None))
            continue
        ext_s, ext_e, _ = hits[key]
        trunc_s, trunc_e = key[1], key[2]
        left_rec = trunc_s - ext_s          # bp regained on the genomic left
        right_rec = ext_e - trunc_e
        left_over = max(0, true_s - ext_s)  # bp beyond the curated end
        right_over = max(0, ext_e - true_e)
        rows.append((seqid, true_s, true_e, orient,
                     (left_rec, right_rec, left_over, right_over)))

    per_member = os.path.join(work, "members.tsv")
    with open(per_member, "w") as out:
        out.write("seqid\ttrue_start\ttrue_end\torient\t"
                  "left_recovered\tright_recovered\tleft_overshoot\tright_overshoot\n")
        for seqid, s, e, o, r in rows:
            if r is None:
                out.write(f"{seqid}\t{s}\t{e}\t{o}\tNA\tNA\tNA\tNA\n")
            else:
                out.write(f"{seqid}\t{s}\t{e}\t{o}\t{r[0]}\t{r[1]}\t{r[2]}\t{r[3]}\n")

    scored = [r[4] for r in rows if r[4] is not None]
    n = len(scored)
    if n:
        # Recovery capped at T per side for the summary fraction; raw
        # overshoot reported separately so over-extension isn't rewarded.
        mean = lambda xs: sum(xs) / len(xs)
        lrec = mean([min(x[0], t) / t for x in scored])
        rrec = mean([min(x[1], t) / t for x in scored])
        lov = mean([x[2] for x in scored])
        rov = mean([x[3] for x in scored])
        full = sum(1 for x in scored
                   if min(x[0], t) / t >= 0.9 and min(x[1], t) / t >= 0.9)
        print(f"{os.path.basename(stk_path)} T={t}: members={n} "
              f"(dropped={len(members) - len(truth)}, missing={len(truth) - n})")
        print(f"  mean recovery  left={lrec:.2f} right={rrec:.2f} "
              f"(fraction of {t}bp; 1.0 = perfect)")
        print(f"  mean overshoot left={lov:.0f}bp right={rov:.0f}bp")
        print(f"  members with >=90% recovery both sides: {full}/{n}")
        print(f"  per-member: {per_member}")
    else:
        print("no members scored (all missing from ram-extend output)")


if __name__ == "__main__":
    main()
