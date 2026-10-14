# SAND-64 — Differential characteristics found by SAT

Reproducible artefacts for the differential cryptanalysis of **SAND-64**.

Main results, all obtained with the SAT model in this repository and re-checked
by an independent transfer-matrix computation:

| Rounds | Minimum weight $-\log_2 p$ | Optimality | Characteristics attaining it |
| --- | --- | --- | --- |
| 12 | **58** ( $p = 2^{-58}$ ) | weight $\le 57$ is UNSAT | **12 288**, exhaustively enumerated and proved complete |
| 13 | **64** ( $p = 2^{-64}$ ) | weight $\le 63$ is UNSAT | one example is published here (not enumerated exhaustively) |

The 12 288 optimal 12-round characteristics involve 352 distinct input
differences and 352 distinct output differences, forming 6 144 distinct pairs
$(\Delta P, \Delta C)$. **Every such pair is realized by exactly two
characteristics**, so the corresponding 12-round differential has probability at
least $2 \cdot 2^{-58} = 2^{-57}$.

Only four round-weight vectors $(w_1,\dots,w_{12})$ are realizable; they form
two time-reversal pairs:

| $(w_1,\dots,w_{12})$ | count |
| --- | --- |
| (4, 2, 0, 2, 4, 9, 9, 8, 6, 4, 6, 4) | 5 376 |
| (4, 6, 4, 6, 8, 9, 9, 4, 2, 0, 2, 4) | 5 376 |
| (2, 7, 4, 7, 8, 9, 9, 4, 2, 0, 2, 4) | 768 |
| (4, 2, 0, 2, 4, 9, 9, 8, 7, 4, 7, 2) | 768 |

## The joint DDT of SSb — why a naive model is wrong

A SAND round applies two branch functions $G_0$ and $G_1$ to the **same** nibble
of the round input. Their output differences must therefore be drawn from the
*joint* differential distribution table of the synthetic S-box

$$\mathrm{SSb}(x) = G_0(x) \,\|\, G_1(x).$$

Taking the product of the two marginal DDTs admits **157** triples
$(\delta_{\mathrm{in}}, \delta_{G_0}, \delta_{G_1})$, whereas only **97** are
jointly realizable. A model built from independent DDTs consequently produces
"characteristics" containing round transitions of probability exactly zero.
For example, the marginal model accepts input nibble difference `0x6` together
with $(\delta_{G_0}, \delta_{G_1}) = (0\mathrm{x}6, 0\mathrm{x}2)$, while the
joint entry is `DDT_SSb[0x6][0x62] = 0`.

The CNF encoding of the joint DDT used here has 58 clauses over 15 variables and
its exactness is *proved*, not sampled: `verify_ssb_model()` enumerates all
$2^{15}$ assignments and compares the models against the table.

## Requirements

* Python 3 (standard library only)
* [Kissat](https://github.com/arminbiere/kissat) on `PATH`, or at `./build/kissat`

Timings below were measured on Ubuntu 24.04.2 LTS, AMD Ryzen 7 5800X, 128 GB RAM.

## Quick verification (recommended for reviewers)

This is the fastest independent check and it does **not** use the SAT solver at
all. For every round transition it counts the exact number of 32-bit round
inputs realizing it by a transfer-matrix product over the eight SSb positions
(so the joint $G_0/G_1$ constraint is enforced on every shared nibble) and
derives the weight as $32 - \log_2(\#)$:

```
python3 verify_characteristics.py
```

It verifies that all 12 288 published 12-round characteristics have exact
weight 58 with no impossible round transition, that they are pairwise distinct,
and that the published 13-round characteristic has exact weight 64. Runtime is
a few minutes.

To check one file on its own:

```
python3 verify_characteristics.py results/trails_w58_all.txt 58
```

## Proving completeness of the 12-round list

`verify_characteristics.py` shows that every published characteristic is valid;
the following shows that **no other one exists**. The 12-round model is built
with "total weight exactly 58" and *no* partitioning constraint, and every
published characteristic is added as a blocking clause:

```
python3 completeness_check.py
```

Expected output: `UNSATISFIABLE` (about 2 minutes). This argument does not
depend on how the enumeration was organised.

## Reproducing the search from scratch

Optimality of the bounds — each command prints `SATISFIABLE` or
`UNSATISFIABLE` together with its runtime:

```
python3 enumerate_w58.py check 58 12     # SATISFIABLE    (~23 s)
python3 enumerate_w58.py check 57 12     # UNSATISFIABLE  (~55 s)  => 58 is minimal
python3 enumerate_w58.py check 64 13     # SATISFIABLE    (~6 s)
python3 enumerate_w58.py check 63 13     # UNSATISFIABLE  (~182 s) => 64 is minimal
```

Full enumeration of the weight-58 12-round characteristics. The search space is
partitioned first by the round-weight vector and then by the set of active SSb
positions in round 1, giving 92 disjoint sub-slices that are solved in parallel;
each sub-slice is enumerated by repeated solving with blocking clauses:

```
python3 enum_by_wvec.py vectors          # phase A: the 4 realizable weight vectors
./run_until_complete.sh 15               # phase B + aggregation + completeness proof
```

`run_until_complete.sh` takes the number of parallel workers as its argument. It
loops until the completeness check returns UNSAT, then compares the regenerated
`trails_w58_all.txt` with the published `results/trails_w58_all.txt`.

Note that the sub-slice partition is only a means of parallelisation: whatever
it produces, the final answer is accepted only once the partition-independent
completeness check passes.

## Key-recovery experiments (Section V-C of the paper)

The directory `experiments/` contains the experiments that validate the
chosen-ciphertext key-recovery attack on 16-round SAND-64. Both scripts are
self-contained (they re-implement SAND-64 and do not import the search code)
and require Python 3 with **NumPy**.

```
cd experiments
python3 cc_check.py      # about 2.5 minutes
python3 cc_random16.py   # about 1.5 minutes
```

`cc_check.py` runs an 8-round version of the attack with independent random
round keys. It keeps the three key-recovery rounds (guessing
`sk^0 & DFDFDFDF` and `sk^1 & 8C008C8C`) and the swap-free last round of the
16-round attack, and replaces the 12-round characteristic by the 4-round
characteristic
`(dx3,dy3) = (00900000,008C0000) -> (dx7,dy7) = (00900000,00010000)`
(weight 8; the last four rounds of the 12-round characteristic read in the
encryption direction, identical to rounds 3-6 of the 16-round attack).
Ciphertext pairs are chosen as in Step 1 of the attack. It checks

* that the chosen ciphertext pairs give the difference `(dx7,dy7)`,
* that every right pair passes the two plaintext-side filters and is counted
  by the right key class,
* that the test depends only on the 37 guessed key bits, and
* the fraction of right pairs counted by wrong key classes close to the right
  one, from which the expected number of wrong classes reaching the threshold
  `tau = 16` and the probability that the right class is not ranked first are
  derived for `lambda = 32` right pairs.

`cc_random16.py` decrypts `2^26` chosen-ciphertext pairs through 16 rounds and
measures the pass rates of the two filters (random-difference heuristic).

The outputs obtained for the paper are stored in
`experiments/cc_check_output.txt` and `experiments/cc_random16_output.txt`;
the scripts use a fixed seed and reproduce them exactly.

## Files

| File | Role |
| --- | --- |
| `SAND-64_20260729.py` | SAT model of SAND-64 (joint SSb DDT, cardinality constraints, Kissat interface). Running it directly performs a single bounded search. |
| `transfer_matrix.py` | Exact transfer-matrix computation of round-transition counts, independent of the SAT model (used by the verification and aggregation scripts). |
| `trail_io.py` | Reader/writer for the characteristic file format. |
| `verify_characteristics.py` | **Standalone verification of the published characteristics.** |
| `completeness_check.py` | **Proof that the 12-round list is complete** (single UNSAT call). |
| `enumerate_w58.py` | Exactly-$k$ cardinality constraint, single-slice enumeration, `check` mode for the bounds. |
| `enum_by_wvec.py` | Phase A (realizable round-weight vectors) and per-vector enumeration. |
| `enum_fine.py` | Phase B: 92 sub-slices enumerated in parallel, resumable. |
| `aggregate_w58.py` | De-duplication, statistics and transfer-matrix verification of the raw enumeration output. |
| `run_until_complete.sh` | Driver: enumerate → aggregate → prove completeness, repeat until UNSAT. |
| `results/trails_w58_all.txt` | All 12 288 optimal 12-round characteristics. |
| `results/characteristic_13r_w64.txt` | An optimal 13-round characteristic. |
| `results/weight_vectors_58.json` | The four realizable round-weight vectors. |
| `results/results_w58_summary.txt` | Summary of the enumeration (in Japanese). |
| `experiments/cc_check.py` | 8-round experiment of the chosen-ciphertext key-recovery attack (Section V-C). |
| `experiments/cc_random16.py` | Filter pass rates on 16 rounds (random-difference heuristic). |
| `experiments/*_output.txt` | Outputs of the two experiments reported in the paper. |

## File format

```
# characteristic 1
  0 0012 0000 0019 0000
  1 0002 0000 0012 0000
  ...
 12 008A 0000 0090 0000
```

Each line is `r  dX^r  dY^r`, where `(dX^0, dY^0)` is the plaintext difference
and `(dX^R, dY^R)` the ciphertext difference. Bit 31 of `dX^r` is the most
significant bit of the left half of the state; nibble `j` of the SSb layer is
formed from bits `31-j`, `23-j`, `15-j`, `7-j`, as in the SAND specification.
In the forward model `dY^r = dX^{r-1}`, so a characteristic is determined by
`dY^0` together with `dX^0, ..., dX^R`.

## Example — an optimal 12-round characteristic

| r | dX^r | dY^r | log2(p) |
| --- | --- | --- | --- |
| 0 | 00120000 | 00190000 | |
| 1 | 00020000 | 00120000 | −4 |
| 2 | 00000000 | 00020000 | −6 |
| 3 | 00020000 | 00000000 | −6 |
| 4 | 00120000 | 00020000 | −8 |
| 5 | 00190200 | 00120000 | −12 |
| 6 | 00C10200 | 00190200 | −21 |
| 7 | 002B0000 | 00C10200 | −30 |
| 8 | 002A0000 | 002B0000 | −38 |
| 9 | 00500000 | 002A0000 | −44 |
| 10 | 00070000 | 00500000 | −48 |
| 11 | 00900000 | 00070000 | −54 |
| 12 | 008A0000 | 00900000 | −58 |
