"""Ejercicio 1: 10 secuencias de una misma familia, su matriz de puntaje /
identidad, si es simetrica, y que cambia segun el puntaje que se use."""

import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import matrices as mx
import plots
from blast import run_blast
from ncbi import fetch_family, load_metadata


def sequences():
    """The records and their ids, in the fixed order every matrix uses."""
    records = fetch_family()
    return records, [r.id for r in records]


# ---------------------------------------------------------------------------
# i) 10 secuencias de una misma familia
# ---------------------------------------------------------------------------

def i():
    records, ids = sequences()
    meta = load_metadata()

    print(f"\n{len(records)} sequences, beta-globin, one per species\n")
    print(f"{'id':<24s}{'acc':<9s}{'entry':<13s}{'len':>5s}  why this entry")
    for r in records:
        m = meta[r.id]
        print(f"{r.id:<24s}{m['accession']:<9s}{m['entry_name']:<13s}"
              f"{m['length']:>5d}  {m['chosen_because']}")

    lengths = [len(r.seq) for r in records]
    print(f"\nlengths: {min(lengths)}-{max(lengths)} aa "
          f"(mean {sum(lengths)/len(lengths):.1f})")
    # Worth a line in the report: the length spread is what makes the raw bit
    # score unusable as a distance without normalising, and what makes
    # BLAST's pident disagree with the coverage-corrected identity.


# ---------------------------------------------------------------------------
# ii) alinearlas y obtener la matriz de puntaje / identidad
# ---------------------------------------------------------------------------

def ii():
    _, ids = sequences()
    hits = run_blast()
    M = mx.build(hits, ids)

    absent = mx.missing_pairs(hits, ids)
    print(f"\n{len(hits)} ordered pairs aligned, {len(absent)} not aligned at all")
    for a, b in absent:
        print(f"  no hit: {a} -> {b}")

    mx.show(M["pident"], ids, "i) " + mx.LABELS["pident"], "{:8.2f}")
    mx.show(M["gident"], ids, "i) " + mx.LABELS["gident"], "{:8.2f}")
    mx.show(M["bitscore"], ids, "ii) " + mx.LABELS["bitscore"], "{:8.1f}")
    mx.show(M["evalue"], ids, "iii) " + mx.LABELS["evalue"], "{:8.1e}")
    return M, ids


# ---------------------------------------------------------------------------
# Analice si la matriz es (o no) simetrica
# ---------------------------------------------------------------------------

def symmetry():
    _, ids = sequences()
    M = mx.build(run_blast(), ids)

    print("\nSymmetry of the raw matrices (off-diagonal pairs, i vs j against j vs i)\n")
    print(f"{'metric':<12s}{'pares distintos':>18s}{'max |dif|':>14s}{'dif relativa media':>22s}")
    for m in mx.METRICS:
        n_diff, max_abs, mean_rel = mx.symmetry(M[m])
        print(f"{m:<12s}{n_diff:>12d} / 45{max_abs:>14.4g}{mean_rel:>21.2%}")

    # The interesting cases printed in full, since "the matrix is roughly
    # symmetric" is not an answer: the consigna asks whether it is, and the
    # honest answer is that it depends on the metric.
    for m in mx.METRICS:
        n_diff, _, _ = mx.symmetry(M[m])
        if n_diff == 0:
            continue
        print(f"\n{m}: the pairs that disagree between the two directions")
        worst = []
        for a in range(len(ids)):
            for b in range(a + 1, len(ids)):
                x, y = M[m][a, b], M[m][b, a]
                if np.isfinite(x) and np.isfinite(y) and x != y:
                    scale = max(abs(x), abs(y)) or 1.0
                    worst.append((abs(x - y) / scale, ids[a], ids[b], x, y))
        for rel, a, b, x, y in sorted(worst, reverse=True)[:5]:
            print(f"  {a:>22s} -> {b:<22s} {x:>12.4g}   "
                  f"reverso {y:<12.4g} ({rel:.1%})")

    print("""
Why: identity and bit score are properties of the alignment, and BLAST's
alignment of (i, j) is essentially the alignment of (j, i) -- what little
asymmetry is left comes from composition-based statistics, which rescale the
substitution matrix using the *query's* amino acid composition, and from the
heuristic seeding, which indexes the subject and scans the query.

The E-value is different in kind. E = K*m*n*exp(-lambda*S): it is proportional
to the query length m and to the database length n. Swapping query and subject
swaps m and n, so unless the two sequences have the same length the E-value
*must* differ, and the difference grows with the length ratio. That is not
noise, it is the definition. An E-value is a statement about a search, not
about a pair, which is the reason it is a bad candidate for a distance.
""")


# ---------------------------------------------------------------------------
# Compare matrices utilizando diferentes puntajes
# ---------------------------------------------------------------------------

def compare():
    _, ids = sequences()
    M = mx.build(run_blast(), ids)

    S = {m: mx.symmetrize(M[m], m) for m in mx.METRICS}
    D = {m: mx.to_distance(S[m], m) for m in mx.METRICS}

    for m in mx.METRICS:
        mx.show(D[m], ids, f"Distancias desde {mx.LABELS[m]}", "{:8.3f}", half=True)

    print("\nCorrelacion de rangos (Spearman) entre las matrices de distancia")
    print("      " + "".join(f"{m[:8]:>10s}" for m in mx.METRICS))
    for a in mx.METRICS:
        row = "".join(f"{mx.spearman(mx.upper(D[a]), mx.upper(D[b])):>10.3f}"
                      for b in mx.METRICS)
        print(f"{a[:6]:<6s}{row}")

    # Spearman is the summary; the nearest-neighbour table is the consequence.
    # The guide tree's first fusion is the global minimum of the matrix, so
    # two metrics can correlate at 0.99 and still start the tree differently.
    print("\nVecino mas cercano de cada secuencia, por metrica")
    nn = {m: mx.nearest_neighbours(D[m], ids) for m in mx.METRICS}
    width = max(len(s) for s in ids)
    print(" " * (width + 1) + "".join(f"{m[:14]:>16s}" for m in mx.METRICS))
    for name in ids:
        cells = "".join(f"{nn[m][name][:14]:>16s}" for m in mx.METRICS)
        flag = "" if len({nn[m][name] for m in mx.METRICS}) == 1 else "   <- disagree"
        print(f"{name:<{width}s} {cells}{flag}")

    for m in mx.METRICS:
        i, j = min(((a, b) for a in range(len(ids)) for b in range(len(ids)) if a < b),
                   key=lambda p: D[m][p])
        print(f"\nprimer par que uniria el arbol guia segun {m}: "
              f"{ids[i]} + {ids[j]}  (d = {D[m][i, j]:.4f})")

    plots.heatmaps(D, ids)
    plots.scatter(D, ids)
    return D, ids
