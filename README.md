# SAND — Differential attacks found by SAT

Reproducible artefacts for

> Nobuyuki Sugio, *Differential Attacks on Lightweight Block Cipher Family SAND*
> (Hokkaido University of Science)

Both members of the family are covered by one SAT model, differing only in the
parameters of Table 3 of the paper. Every headline number below is reproduced by
the search code in this repository **and** re-checked by an independent
re-implementation that shares no code with the search.

| | Distinguisher | Key recovery | Data | Time | Memory |
| --- | --- | --- | --- | --- | --- |
| **SAND-64** | 12 rounds, $p = 2^{-57}$ | **17 rounds** | $2^{63}$ | $2^{83.5}$ | $2^{58}$ |
| **SAND-128** | 19 rounds, $p = 2^{-120}$ | **23 rounds** | $2^{126}$ | $2^{122.5}$ | $2^{51}$ |

The 12-round SAND-64 distinguisher reaches $2^{-57}$ rather than $2^{-58}$
because the 12 288 optimal characteristics of weight 58 collapse onto 6 144
distinct $(\Delta P, \Delta C)$ pairs, **each realized by exactly two
characteristics**.

## The joint DDT of SSb — why a naive model is wrong

A SAND round applies two branch functions $G_0$ and $G_1$ to the **same** nibble
of the round input. Their output differences must therefore be drawn from the
*joint* differential distribution table of the synthetic S-box

$$\mathrm{SSb}(x) = G_0(x) \,\|\, G_1(x).$$

Taking the product of the two marginal DDTs admits **157** triples
$(\delta_{\mathrm{in}}, \delta_{G_0}, \delta_{G_1})$, whereas only **97** are
jointly realizable. A model built from independent DDTs consequently produces
"characteristics" containing round transitions of probability exactly zero — for
example input nibble difference `0x6` with
$(\delta_{G_0}, \delta_{G_1}) = (0\mathrm{x}6, 0\mathrm{x}2)$, while the joint
entry is `DDT_SSb[0x6][0x62] = 0`. This holds for both variants, since SSb is
the same 4-bit map in each.

The CNF encoding of the joint DDT is *proved* exact, not sampled:
`verify_ssb_model()` enumerates all $2^{15}$ assignments and compares the models
against the table.

## Layout

```
sand64/    SAND-64:  SAT model, enumeration, 17-round key recovery
sand128/   SAND-128: SAT model, enumeration, 23-round key recovery
build/     Kissat binary (not tracked; see Requirements)
```

Each variant directory is self-contained and keeps its own `results/`. All
scripts locate their data relative to their own file, so they can be run from
any working directory.

## Requirements

* Python 3 (standard library only)
* [Kissat](https://github.com/arminbiere/kissat) at `build/kissat` in the
  repository root, or anywhere on `PATH`

Timings below were measured on Ubuntu 24.04.2 LTS, AMD Ryzen 7 5800X, 128 GB RAM.

## Quick verification (recommended for reviewers)

These two scripts are the fastest independent check and they do **not** use the
SAT solver at all. Each re-implements its cipher from the specification in
Section IV of the paper and recomputes every claim — round-transition
probabilities by an exact transfer-matrix product over the SSb positions, the
key-free filters and guessed-key masks by exhaustive bit tracking, and the
counting phase by simulation:

```
python3 sand64/verify_sand64_key_recovery.py     # 17-round attack, 44 checks
python3 sand128/verify_sand128_key_recovery.py   # 23-round attack
```

Both must end with `ALL CHECKS PASSED`. The author's own search modules are
imported only in the final section, to cross-check that the two independent
implementations agree.

## Regression test

CNF generation is fully deterministic, so the SHA-256 of a generated CNF is a
fingerprint of the model. `tests/cnf_fingerprint.py` builds sixteen instances
across both variants — Matsui pruning on and off, the exactly-$k$ enumeration
model, the weight-vector partition, and the key-recovery boundary constraints —
and compares their fingerprints against `tests/cnf_baseline.json`:

```
python3 tests/cnf_fingerprint.py            # ~2 s, no SAT solving
python3 tests/cnf_fingerprint.py --update   # re-record, only after an intended change
```

This is what makes restructuring the code safe: if the fingerprints are
unchanged, the search behaves identically, without re-running any search. It is
sensitive to changes no size check would catch — the plaintext- and
ciphertext-side boundary constraints produce CNFs with identical variable and
clause counts but different fingerprints, so swapping them fails the test.

## Paper ↔ command

| Result in the paper | Command | Data file |
| --- | --- | --- |
| Minimum weights $B[r]$, SAND-64 | `python3 sand64/matsui_sand64.py` | `sand64/bounds_sand64.json` |
| Minimum weights $B[r]$, SAND-128 | `python3 sand128/matsui_sand128.py` | `sand128/bounds_sand128.json` |
| 12-round characteristic, SAND-64 | `python3 sand64/enumerate_w58.py check 58 12` | `sand64/results/trails_w58_all.txt` |
| Optimality of 58 | `python3 sand64/enumerate_w58.py check 57 12` → UNSAT | |
| Round-weight vectors of the 12 288 | `python3 sand64/enum_by_wvec.py vectors` | `sand64/results/weight_vectors_58.json` |
| Completeness of the 12-round list | `python3 sand64/completeness_check.py` → UNSAT | |
| 19-round characteristic, SAND-128 | `python3 sand128/enumerate_sand128.py 19 120` | `sand128/results/trails_sand128_r19_w120.txt` |
| Key-free filters (Eq. 9–12) | `python3 sand128/filters_sand128.py validate` | |
| Guessed round-key bits $k$ | `python3 sand128/keybits_sand128.py` | |
| **All numbers of the 23-round attack** | `python3 sand128/attack_numbers.py` | |
| Threshold $\tau$ and success probability | `python3 sand128/success_probability_sand128.py` | |
| Maximum attackable rounds (Discussion) | `python3 sand128/maximize_rounds.py` | |

## SAND-64: enumeration and its completeness

`verify_characteristics.py` shows that every published characteristic is valid;
`completeness_check.py` shows that **no other one exists**. The 12-round model is
rebuilt with "total weight exactly 58" and *no* partitioning constraint, and
every published characteristic is added as a blocking clause:

```
python3 sand64/verify_characteristics.py     # a few minutes
python3 sand64/completeness_check.py         # expects UNSATISFIABLE, ~2 min
```

Expected `UNSATISFIABLE` — an argument that does not depend on how the
enumeration was organised.

To reproduce the enumeration itself, the search space is partitioned first by
round-weight vector and then by the set of active SSb positions in round 1,
giving 92 disjoint sub-slices solved in parallel:

```
python3 sand64/enum_by_wvec.py vectors       # phase A: the 4 realizable vectors
cd sand64 && ./run_until_complete.sh 15      # phase B + aggregation + completeness
```

`run_until_complete.sh` takes the number of parallel workers as its argument and
loops until the completeness check returns UNSAT, then compares the regenerated
list with the published `results/trails_w58_all.txt`.

Results, for reference:

| $(w_1,\dots,w_{12})$ | count |
| --- | --- |
| (4, 2, 0, 2, 4, 9, 9, 8, 6, 4, 6, 4) | 5 376 |
| (4, 6, 4, 6, 8, 9, 9, 4, 2, 0, 2, 4) | 5 376 |
| (2, 7, 4, 7, 8, 9, 9, 4, 2, 0, 2, 4) | 768 |
| (4, 2, 0, 2, 4, 9, 9, 8, 7, 4, 7, 2) | 768 |

## SAND-128: enumeration

Same strategy, ported to the 16-column state:

```
python3 sand128/enum_by_wvec_sand128.py vectors [workers]   # phase A
python3 sand128/enum_by_wvec_sand128.py run     [workers]   # phase B
```

Unlike SAND-64, the 19-round weight-120 list is **not** claimed to be
exhaustive; 151 distinct characteristics of weight 120 were found and every one
of them was re-checked by the transfer-matrix computation.

## Files

### `sand64/`

| File | Role |
| --- | --- |
| `SAND-64_20260729.py` | SAT model (joint SSb DDT, cardinality constraints, Kissat interface) |
| `matsui_sand64.py` | Matsui-style lower-bound pruning; computes $B[r]$ incrementally |
| `enumerate_w58.py` | Exactly-$k$ cardinality constraint, single-slice enumeration, `check` mode |
| `enum_by_wvec.py` | Phase A (realizable round-weight vectors) and per-vector enumeration |
| `enum_fine.py` | Phase B: 92 sub-slices enumerated in parallel, resumable |
| `aggregate_w58.py` | De-duplication, statistics, transfer-matrix verification |
| `run_until_complete.sh` | Driver: enumerate → aggregate → prove completeness, until UNSAT |
| `transfer_matrix_sand64.py` | Exact round-transition counts by a transfer-matrix product, independent of the SAT model; used by the enumeration and verification scripts |
| `verify_characteristics.py` | Standalone verification of the published characteristics |
| `completeness_check.py` | Proof that the 12-round list is complete (single UNSAT call) |
| `verify_sand64_key_recovery.py` | **Independent verification of the 17-round attack** |
| `trail_io.py` | Reader/writer for the characteristic file format |
| `results/` | The 12 288 characteristics, the 13-round example, the weight vectors |

### `sand128/`

| File | Role |
| --- | --- |
| `SAND-128_..._20260805.py` | SAT model |
| `matsui_sand128.py` | Matsui-style lower-bound pruning; computes $B[r]$ |
| `run_sand128_bound.py` | Driver for a single bounded search |
| `enumerate_sand128.py` | All-SAT enumeration for given rounds and weight bound |
| `enum_by_wvec_sand128.py` | Enumeration partitioned by round-weight vector |
| `keyrecovery_sand128.py` | Characteristic search under the key-recovery constraints |
| `filters_sand128.py` | Key-free filter conditions by truncated-difference propagation |
| `keybits_sand128.py` | Exact count of the guessed round-key bits $k$ |
| `attack_numbers.py` | **Derives every number of the 23-round attack** from the characteristic data and checks it against the paper |
| `maximize_rounds.py` | Maximizes $r_{in} + R + r_{out}$ under the key-size constraint |
| `success_probability_sand128.py` | Threshold $\tau$ and success probability of the counting phase |
| `verify_sand128_key_recovery.py` | **Independent verification of the 23-round attack** |
| `results/` | The 19-round weight-120 characteristics |

## Characteristic file format

```
# characteristic 1
  0 0012 0000 0019 0000
  1 0002 0000 0012 0000
  ...
 12 008A 0000 0090 0000
```

Each line is `r  dX^r  dY^r`, where $(\Delta X^0, \Delta Y^0)$ is the plaintext
difference and $(\Delta X^R, \Delta Y^R)$ the ciphertext difference. Bit 31
(resp. 63) of $\Delta X^r$ is the most significant bit of the left half of the
state; nibble $j$ of the SSb layer is formed from bits `31-j`, `23-j`, `15-j`,
`7-j` for SAND-64 and `63-j`, `47-j`, `31-j`, `15-j` for SAND-128, as in the
specification. In the forward model $\Delta Y^r = \Delta X^{r-1}$, so a
characteristic is determined by $\Delta Y^0$ together with
$\Delta X^0, \dots, \Delta X^R$.

## Branches and tags

`main` holds the journal version, covering both variants. Conference versions of
this work are kept out of the way:

| Ref | Contents |
| --- | --- |
| `main` | Journal version: SAND-64 (17 rounds) and SAND-128 (23 rounds) |
| `paper/candarw2026` | CANDAR-W 2026 version: SAND-64 only, 16-round attack |
| `candarw2026-submitted` | The exact state submitted to CANDAR-W 2026 |

## License

MIT — see `LICENSE`.
