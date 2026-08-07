#!/bin/bash
# Re-run the whole enumeration of the weight-58 12-round characteristics from
# scratch and iterate until the completeness check reports UNSATISFIABLE.
#
# Each iteration
#   1. clears the per-sub-slice markers and re-enumerates all 92 sub-slices;
#      characteristics already known are supplied as blocking clauses, but only
#      those that can belong to the sub-slice in question,
#   2. copies the results into trails_parts/ AFTER the worker pool has finished
#      (copying while the workers are still writing loses characteristics),
#   3. aggregates and de-duplicates everything into trails_w58_all.txt,
#   4. runs the partition-independent completeness check on that file.
#
# The loop stops as soon as the check is UNSAT, which proves the enumeration
# complete.  A run that starts from an empty trails_parts/ needs several
# iterations; starting from the published results it converges immediately.
set -u
cd "$(dirname "$0")"
mkdir -p trails_parts trails_w58_fine

for iter in $(seq 1 8); do
    echo "===================== iteration $iter ====================="
    rm -f trails_w58_fine/*.done trails_w58_fine/*.txt
    python3 -u enum_fine.py "${1:-15}" > "phase_iter${iter}.log" 2>&1
    echo "--- enumeration finished"
    tail -4 "phase_iter${iter}.log"

    # Snapshot only after the pool has terminated.
    for f in trails_w58_fine/*.txt; do
        [ -f "$f" ] && cp "$f" "trails_parts/i${iter}_$(basename "$f")"
    done

    python3 -u aggregate_w58.py > "aggregate_iter${iter}.log" 2>&1
    grep -E "total distinct|written to" "aggregate_iter${iter}.log"

    python3 -u completeness_check.py trails_w58_all.txt > "completeness_iter${iter}.log" 2>&1
    echo "--- completeness check:"
    tail -3 "completeness_iter${iter}.log"

    if grep -q "^UNSATISFIABLE" "completeness_iter${iter}.log"; then
        echo "=== enumeration is COMPLETE after $iter iteration(s) ==="
        echo "=== compare with the published list:"
        diff <(grep -v '^#' trails_w58_all.txt) \
             <(grep -v '^#' results/trails_w58_all.txt) > /dev/null \
             && echo "    identical to results/trails_w58_all.txt" \
             || echo "    DIFFERS from results/trails_w58_all.txt (compare as sets, order may vary)"
        exit 0
    fi
done
echo "=== still incomplete after all iterations ==="
exit 1
