"""Ejercicio 3) Construya un algoritmo que dados los AP enre las secuencias del objetivo 1 y el arbol
guía genere el AM y lo muestre en pantalla. Compare sus resultados con los obtenidos utilizando
algún programa on-line.
Ayuda: para ir ubicando las bases de cada secuencia en el AM utilice un indice “unico” para el AM"""

import os
import shutil
import subprocess

import matrices as mx
import msa
import plots
from exercises.ex1 import family
from sequences import DATA_DIR, run_blast

# The guide tree of the MSA: the global identity one, the one exercise 2 uses
# as reference in all its figures.
GUIDE = "gident"
MSA_FASTA = os.path.join(DATA_DIR, "msa.fasta")

# The programs it is compared against. MAFFT rather than Clustal Omega because
# `brew install mafft` installs it (Clustal Omega is not in Homebrew). Two
# configurations: FFT-NS-2 is progressive, like our method (guide tree and
# merges), and L-INS-i is iterative and the most accurate it offers.
REFERENCES = {
    "mafft_fft2":  ["--retree", "2", "--maxiterate", "0"],
    "mafft_linsi": ["--localpair", "--maxiterate", "1000"],
}


def inputs():
    """The sequences, their ids, BLAST's APs and the raw matrices."""
    records, ids = family()
    seqs = {r.id: str(r.seq) for r in records}
    hits = run_blast()
    return seqs, ids, hits, mx.build(hits, ids)


def tree_of(raw, ids, metric):
    D = mx.to_distance(mx.symmetrize(raw[metric], metric), metric)
    return mx.build_tree(D, ids)


def reference_path(name):
    return os.path.join(DATA_DIR, f"ref_{name}.fasta")


def references(force=False):
    """{name: MSA} of each external program, cached in data/ref_*.fasta.

    Any `data/ref_<name>.fasta` added by hand (for example Clustal Omega's
    output, as an aligned FASTA) is compared like the others."""
    fasta = os.path.join(DATA_DIR, "globins.fasta")
    for name, args in REFERENCES.items():
        path = reference_path(name)
        if os.path.exists(path) and not force:
            continue
        if shutil.which("mafft") is None:
            print(f"  (sin mafft instalado: se saltea {name}; `brew install mafft`)")
            continue
        cmd = ["mafft", "--quiet", *args, fasta]
        print("$ " + " ".join(cmd))
        out = subprocess.run(cmd, capture_output=True, text=True, check=True).stdout
        with open(path, "w") as f:
            f.write(out)

    found = {}
    for file in sorted(os.listdir(DATA_DIR)):
        if file.startswith("ref_") and file.endswith(".fasta"):
            found[file[4:-6]] = msa.read_fasta(os.path.join(DATA_DIR, file))
    return found


# ---------------------------------------------------------------------------
# a) The MSA, shown on screen
# ---------------------------------------------------------------------------

def a():
    """Construya un algoritmo que dados los AP enre las secuencias del objetivo 1 y el arbol
    guía genere el AM y lo muestre en pantalla."""
    seqs, ids, hits, raw = inputs()
    table = msa.pairwise_table(hits, seqs)
    root = tree_of(raw, ids, GUIDE)

    print(f"\nArbol guia ({GUIDE}):\n{root.newick};")

    alignment, log = msa.build(root.nested, seqs, table)

    print("\nFusiones, en el orden en que las marca el arbol")
    print(f"{'#':>2s}  {'izquierda':>9s} + {'derecha':<9s}  {'par ancla (mejor AP entre los dos lados)':<52s}{'bit score':>9s}")
    for k, step in enumerate(log, 1):
        x, y = step["anchor"]
        print(f"{k:>2d}  {len(step['left']):>6d} seq + {len(step['right']):<3d} seq  "
              f"{x + ' ~ ' + y:<52s}{step['score']:>9.0f}")

    print(f"\nAlineamiento multiple ({len(next(iter(alignment.values())))} columnas)\n")
    msa.show(alignment)
    print("(* columna identica, : residuos muy similares, . algo similares)")

    st = msa.stats(alignment)
    print(f"\n{st['length']} columnas, {st['identical']} identicas, "
          f"{st['gaps']} caracteres de gap, identidad media por par {st['mean_identity']:.1%}")

    kept, per_pair = msa.consistency(alignment, table)
    worst = sorted(per_pair.items(), key=lambda kv: kv[1])[:3]
    print(f"consistencia con los 45 AP: {kept:.2%} de los pares de residuos de los AP "
          f"estan en el AM; los peores: "
          + ", ".join(f"{x} ~ {y} ({v:.1%})" for (x, y), v in worst))

    msa.write_fasta(alignment, MSA_FASTA)
    print(f"guardado en {MSA_FASTA}")
    plots.msa(alignment, f"Multiple alignment from the guide tree ({GUIDE})",
              symbols=msa.conservation(alignment))
    return alignment


# ---------------------------------------------------------------------------
# b) Comparison with an external program
# ---------------------------------------------------------------------------

def _ranges(positions):
    """{3, 4, 5, 9} -> '3-5, 9' (posiciones 1-based)."""
    out, run = [], []
    for p in sorted(positions):
        if run and p == run[-1] + 1:
            run.append(p)
        else:
            if run:
                out.append(run)
            run = [p]
    if run:
        out.append(run)
    return ", ".join(f"{r[0]}-{r[-1]}" if len(r) > 1 else f"{r[0]}" for r in out)


def b():
    """Compare sus resultados con los obtenidos utilizando algún programa on-line."""
    seqs, ids, hits, raw = inputs()
    mine, _ = msa.build(tree_of(raw, ids, GUIDE).nested, seqs, msa.pairwise_table(hits, seqs))

    print("\nProgramas externos")
    refs = references()
    if not refs:
        print("no hay ningun data/ref_*.fasta ni mafft para generarlo")
        return

    every = {"arbol guia + AP de BLAST (el nuestro)": mine,
             **{name: aln for name, aln in refs.items()}}

    print(f"\n{'':<42s}{'cols':>6s}{'identicas':>11s}{'gaps':>6s}{'id. media':>11s}{'SP (BLOSUM62)':>15s}")
    for name, aln in every.items():
        st = msa.stats(aln)
        print(f"{name:<42s}{st['length']:>6d}{st['identical']:>11d}{st['gaps']:>6d}"
              f"{st['mean_identity']:>11.1%}{msa.sum_of_pairs(aln):>15.0f}")

    print("\nPares de residuos alineados en comun con cada referencia")
    print(f"{'referencia':<20s}{'de la ref. que tenemos':>26s}{'de lo nuestro que esta en la ref.':>36s}")
    for name, aln in refs.items():
        recall, precision = msa.agreement(mine, aln)
        n = len(msa.aligned_pairs(aln))
        print(f"{name:<20s}{recall:>25.2%}{precision:>36.2%}   ({n} pares en la ref.)")

    names = list(refs)
    if len(names) > 1:
        recall, precision = msa.agreement(refs[names[0]], refs[names[1]])
        print(f"\n(entre las dos referencias: {recall:.2%} / {precision:.2%})")

    # Where they disagree. A pair that is in one and not in the other involves
    # two residues; they are listed per sequence, by position in the sequence
    # (not in the column, which changes from one MSA to another).
    for name, aln in refs.items():
        diff = msa.aligned_pairs(mine) ^ msa.aligned_pairs(aln)
        involved = {n: set() for n in seqs}
        for (x, i, y, j) in diff:
            involved[x].add(i + 1)
            involved[y].add(j + 1)
        print(f"\nResiduos en pares que no coinciden con {name} ({len(diff)} pares)")
        for n in ids:
            if involved[n]:
                print(f"  {n:<24s}{_ranges(involved[n])}")
