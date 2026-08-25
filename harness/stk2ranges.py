#!/usr/bin/env python3
"""Convert a RepeatModeler seed alignment (Stockholm) to a RAMExtend ranges TSV.

Mirrors the logic of RepeatAfterMe util/extend-stk.pl:
  - sequence names of the form  <seqid>:<start>-<end>_<orient>  (1-based, inclusive)
  - start is converted to 0-based half-open
  - left/right extendable flags derived from leading/trailing '.' runs:
    up to 10bp of padding from the alignment edge still counts as extendable
  - sequences extendable on neither side are dropped

Usage: stk2ranges.py input.stk output.tsv
"""
import re
import sys


def main(stk_path, tsv_path):
    n_out = 0
    with open(stk_path) as fh, open(tsv_path, "w") as out:
        for line in fh:
            line = line.rstrip("\n")
            if not line or line.startswith("#") or line.startswith("//"):
                continue
            name, seq = line.split(None, 1)
            seq = seq.strip()
            m = re.match(r"^(\S+):(\d+)-(\d+)_([+-])$", name)
            if not m:
                sys.exit(f"Unparseable sequence name: {name}")
            seqid, start, end, orient = m.group(1), int(m.group(2)), int(m.group(3)), m.group(4)
            start -= 1  # to 0-based half-open
            left_ext = 1 if re.match(r"^[.]{0,10}[^.]", seq) else 0
            right_ext = 1 if re.search(r"[^.][.]{0,10}$", seq) else 0
            if not (left_ext or right_ext):
                continue
            out.write(f"{seqid}\t{start}\t{end}\t{left_ext}\t{right_ext}\t{orient}\n")
            n_out += 1
    print(f"{tsv_path}: {n_out} ranges", file=sys.stderr)


if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    main(sys.argv[1], sys.argv[2])
