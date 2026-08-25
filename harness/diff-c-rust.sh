#!/bin/sh
# Differential harness: run the C RAMExtend and the Rust ram-extend on the
# same inputs and compare the -outtsv / -cons / -outfa files byte-for-byte,
# plus the "Extended right/left" stdout lines that Refiner scrapes.
#
# Usage: diff-c-rust.sh [workdir]
#
# Inputs come from c/test/ (ce10 families), and the C binary is built from
# c/ by default, so the harness needs nothing outside this checkout.
set -u

C_BIN=${C_BIN:-$(dirname "$0")/../c/RAMExtend}
RUST_BIN=${RUST_BIN:-$(dirname "$0")/../target/release/ram-extend}

# `cargo test` does not relink the release binary — build explicitly so the
# harness never validates a stale ram-extend.
if [ -z "${SKIP_BUILD:-}" ]; then
    (cd "$(dirname "$0")/.." && cargo build --release) || exit 1
    [ -x "$C_BIN" ] || make -C "$(dirname "$0")/../c" RAMExtend || exit 1
fi
TESTDIR=${TESTDIR:-$(dirname "$0")/../c/test}
WORK=${1:-$(mktemp -d)}
mkdir -p "$WORK"
HARNESS_DIR=$(dirname "$0")

# Parameter sets: mirror extend-stk.pl (L=20000) on the primary run and add a
# default-parameter run. Format: tag:extra-args
PARAM_SETS="extstk:-L 20000 -bandwidth 14 -matrix 20p43g -minimprovement 27
defaults:
refiner:-bandwidth 40 -matrix 25p43g -minimprovement 30
flank:-L 5000 -addflanking 50"

fails=0
total=0

for fam in ce10-fam1 ce10-fam2 ce10-fam3; do
    stk=$TESTDIR/$fam.stk
    tsv=$WORK/$fam.ranges.tsv
    python3 "$HARNESS_DIR/stk2ranges.py" "$stk" "$tsv" 2>/dev/null || {
        echo "FAIL  $fam: stk2ranges conversion failed"; fails=$((fails+1)); continue;
    }
    echo "$PARAM_SETS" | while IFS=: read -r tag extra; do
        :
    done
    # POSIX sh subshell-safe loop over param sets
    for tag in extstk defaults refiner flank; do
        case $tag in
            extstk)   extra="-L 20000 -bandwidth 14 -matrix 20p43g -minimprovement 27" ;;
            defaults) extra="" ;;
            refiner)  extra="-bandwidth 40 -matrix 25p43g -minimprovement 30" ;;
            flank)    extra="-L 5000 -addflanking 50" ;;
        esac
        total=$((total+1))
        base=$WORK/$fam.$tag
        for side in c rust; do
            # The Rust side runs with -ccompat so the three fixed C bugs are
            # reproduced byte-for-byte for this comparison.
            if [ $side = c ]; then bin=$C_BIN; compat=""; else bin=$RUST_BIN; compat="-ccompat"; fi
            # shellcheck disable=SC2086
            "$bin" -twobit "$TESTDIR/ce10.2bit" -ranges "$tsv" \
                -outtsv "$base.$side.tsv" -cons "$base.$side.cons.fa" \
                -outfa "$base.$side.fa" $compat $extra > "$base.$side.log" 2>&1
        done
        ok=1
        for ext in tsv cons.fa fa; do
            if ! cmp -s "$base.c.$ext" "$base.rust.$ext"; then
                echo "FAIL  $fam/$tag: $ext differs ($base.{c,rust}.$ext)"
                ok=0
            fi
        done
        for pat in "Extended right:" "Extended left :"; do
            cline=$(grep "$pat" "$base.c.log")
            rline=$(grep "$pat" "$base.rust.log")
            if [ "$cline" != "$rline" ]; then
                echo "FAIL  $fam/$tag: stdout '$pat' differs: C='$cline' Rust='$rline'"
                ok=0
            fi
        done
        if [ $ok = 1 ]; then
            echo "OK    $fam/$tag"
        else
            fails=$((fails+1))
        fi
    done
done

echo "----"
echo "workdir: $WORK"
if [ $fails = 0 ]; then
    echo "ALL COMPARISONS PASSED"
else
    echo "$fails parameter-set comparisons FAILED"
    exit 1
fi
