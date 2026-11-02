"""Ejercicio 1): utilice lo aprendido en los trabajos prácticos previos (uso de biopython para obtener
secuencias, BLAST y el “parser del mismo”) para:
i) obtener 10 secuencias pertenecientes a una misma familia y
ii) Alinearlas y obtener la correspondiente matriz de puntaje / identidad
Analice si la matriz es (o no) simetrica
Compare matrices utilizando diferentes puntajes que le otorga BLAST (Bit Score, Identidad, E-Value)."""

import numpy as np

import matrices as mx
import plots
from sequences import fetch_family, gident, load_metadata, run_blast


def family():
    """The records and their ids, in the fixed order every matrix uses."""
    records = fetch_family()
    return records, [r.id for r in records]


# ---------------------------------------------------------------------------
# i) 10 secuencias de una misma familia
# ---------------------------------------------------------------------------

def i():
    """i) obtener 10 secuencias pertenecientes a una misma familia"""
    records, ids = family()
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
    """ii) Alinearlas y obtener la correspondiente matriz de puntaje / identidad"""
    _, ids = family()
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
    """Analice si la matriz es (o no) simetrica"""
    _, ids = family()
    M = mx.build(run_blast(), ids)

    print("\nSymmetry of the raw matrices (off-diagonal pairs, i vs j against j vs i)\n")
    print(f"{'metrica':<12s}{'pares distintos':>18s}{'dif. maxima':>14s}{'dif. relativa media':>22s}")
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
                  f"reverse {y:<12.4g} ({rel:.1%})")

    print("""
Why: identity and bit score are properties of the alignment. BLAST's
alignment of (i, j) is the alignment of (j, i) for every pair (qseq/sseq are
identical in both directions), so identity is exactly symmetric. The bit score
is not: the same alignment is scored differently when the query changes,
because blastp rescales the substitution matrix with the *query's* amino acid
composition (-comp_based_stats 2 by default; with -comp_based_stats 0 the bit
score is symmetric too).

The E-value is different in kind. E = K*m*n*exp(-lambda*S): it is proportional
to the query length m and to the database length n. Swapping query and subject
swaps m and n, so unless the two sequences have the same length the E-value
differs; with lengths of 145-148 aa that factor is tiny, and what moves it by
orders of magnitude is the difference in bit score. An E-value is a statement
about a search, not about a pair, which is the reason it is a bad candidate
for a distance.
""")


# ---------------------------------------------------------------------------
# Compare matrices utilizando diferentes puntajes
# ---------------------------------------------------------------------------

def compare():
    """Compare matrices utilizando diferentes puntajes que le otorga BLAST (Bit Score, Identidad, E-Value)."""
    _, ids = family()
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

    plots.scatter(D, ids)
    return D, ids


# ---------------------------------------------------------------------------
# Note: how much of the matrix is the sequences' and how much the scoring's
# ---------------------------------------------------------------------------

MATRICES = [None, "BLOSUM45", "BLOSUM80"]      # None is blastp's default: BLOSUM62


def matrix_sensitivity():
    """Not part of the assignment: how much of the matrices comes from the
    sequences and how much from the substitution matrix (BLOSUM45, 62 and 80)."""
    _, ids = family()
    runs = {}
    for name in MATRICES:
        runs[name or "BLOSUM62"] = run_blast(matrix=name)

    def off_diagonal(hits):
        return [h for h in hits if h["qseqid"] != h["sseqid"]]

    print(f"\n{'matriz':<10s}{'length medio':>14s}{'min':>6s}{'max':>6s}{'ident. global media':>22s}{'cubren todo':>13s}")
    for name, hits in runs.items():
        rows = off_diagonal(hits)
        lengths = [h["length"] for h in rows]
        covers = sum(1 for h in rows if h["qend"] - h["qstart"] + 1 == h["qlen"]
                     and h["send"] - h["sstart"] + 1 == h["slen"])
        print(f"{name:<10s}{sum(lengths)/len(lengths):>14.1f}{min(lengths):>6d}{max(lengths):>6d}"
              f"{sum(gident(h) for h in rows)/len(rows):>21.2f}%{covers:>9d}/{len(rows)}")

    pairs = [("Homo_sapiens", "Danio_rerio"), ("Homo_sapiens", "Xenopus_laevis"),
             ("Gallus_gallus", "Xenopus_laevis"), ("Homo_sapiens", "Bos_taurus"),
             ("Homo_sapiens", "Macaca_mulatta")]
    print(f"\nlength del HSP (nident) para algunos pares, query -> subject")
    print(f"{'par':<34s}" + "".join(f"{name:>16s}" for name in runs))
    for a, b in pairs:
        cells = []
        for hits in runs.values():
            h = next(h for h in hits if h["qseqid"] == a and h["sseqid"] == b)
            cells.append(f"{h['length']} ({h['nident']})")
        print(f"{a + ' -> ' + b:<34s}" + "".join(f"{c:>16s}" for c in cells))

    # the guide tree (gident) under each matrix, against the BLOSUM62 one
    trees = {}
    for name, hits in runs.items():
        M = mx.build(hits, ids)
        D = mx.to_distance(mx.symmetrize(M["gident"], "gident"), "gident")
        trees[name] = mx.build_tree(D, ids)
    base = mx.clades(trees["BLOSUM62"].nested)
    print("\nArbol guia (gident) segun la matriz, contra el de BLOSUM62")
    for name, root in trees.items():
        clades = mx.clades(root.nested)
        only_here = sorted(clades - base, key=len)
        only_base = sorted(base - clades, key=len)
        print(f"  {name:<10s}{len(clades & base)} de {len(base)} clados en comun")
        for c in only_here:
            print(f"      solo con {name}: {{{', '.join(sorted(c))}}}")
        for c in only_base:
            print(f"      solo con BLOSUM62: {{{', '.join(sorted(c))}}}")
