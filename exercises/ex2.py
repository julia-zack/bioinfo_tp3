"""Ejercicio 2): Programe un algoritmo que dada una matriz de distancias construya un árbol guía.
Utilice el mismo para trabajar con las diferentes matrices obtenidas a partir de los puntajes de
BLAST y analice los resultados de manera comparativa."""

import matrices as mx
import plots
from sequences import run_blast
from exercises.ex1 import family

def a():
    """Programe un algoritmo que dada una matriz de distancias construya un árbol guía.
    Utilice el mismo para trabajar con las diferentes matrices obtenidas a partir de los puntajes de
    BLAST"""
    # Builds and draws the guide tree of each metric.
    _, ids = family()
    raw = mx.build(run_blast(), ids)

    trees = {}
    for metric in mx.METRICS:
        D = mx.to_distance(mx.symmetrize(raw[metric], metric), metric)
        newick = mx.newick(D, ids)
        trees[metric] = newick
        print(f"\n{metric}\n{newick}")

        plots.guide_tree(newick, D, ids,
                         f"Guide tree -- {mx.NAMES[metric]}",
                         f"ex2_tree_{metric}.png")
    return trees


# The species tree these ten sequences should recover (the one in the
# docstring of sequences.py): fish outside, then frog, then chicken, then
# mammals.
EXPECTED_NEWICK = ("(((((((Homo_sapiens,Pan_troglodytes),Macaca_mulatta),"
                   "((Mus_musculus,Rattus_norvegicus),Oryctolagus_cuniculus)),Bos_taurus),"
                   "Gallus_gallus),Xenopus_laevis),Danio_rerio);")

# Its clades, named.
EXPECTED = {
    "primates (humano, chimpance, macaco)":
        {"Homo_sapiens", "Pan_troglodytes", "Macaca_mulatta"},
    "humano + chimpance": {"Homo_sapiens", "Pan_troglodytes"},
    "roedores (raton, rata)": {"Mus_musculus", "Rattus_norvegicus"},
    "roedores + conejo": {"Mus_musculus", "Rattus_norvegicus", "Oryctolagus_cuniculus"},
    "primates + roedores + conejo": {"Homo_sapiens", "Pan_troglodytes", "Macaca_mulatta",
                                      "Mus_musculus", "Rattus_norvegicus", "Oryctolagus_cuniculus"},
    "mamiferos (con la vaca)": {"Homo_sapiens", "Pan_troglodytes", "Macaca_mulatta",
                                 "Mus_musculus", "Rattus_norvegicus", "Oryctolagus_cuniculus",
                                 "Bos_taurus"},
    "amniotas (mamiferos + pollo)": {"Homo_sapiens", "Pan_troglodytes", "Macaca_mulatta",
                                      "Mus_musculus", "Rattus_norvegicus", "Oryctolagus_cuniculus",
                                      "Bos_taurus", "Gallus_gallus"},
    "tetrapodos (todo menos el pez)": {"Homo_sapiens", "Pan_troglodytes", "Macaca_mulatta",
                                        "Mus_musculus", "Rattus_norvegicus", "Oryctolagus_cuniculus",
                                        "Bos_taurus", "Gallus_gallus", "Xenopus_laevis"},
}


def b():
    """analice los resultados de manera comparativa."""
    # Compares the four topologies with each other and with the species tree.
    _, ids = family()
    raw = mx.build(run_blast(), ids)

    clades, newicks = {}, {}
    for metric in mx.METRICS:
        D = mx.to_distance(mx.symmetrize(raw[metric], metric), metric)
        clades[metric] = mx.clades(mx.build_tree(D, ids).nested)
        newicks[metric] = mx.newick(D, ids)

    print("\nClados del arbol de especies que recupera cada metrica")
    print(f"{'':<34s}" + "".join(f"{m:>10s}" for m in mx.METRICS))
    for label, leaves in EXPECTED.items():
        cells = "".join(f"{'si' if frozenset(leaves) in clades[m] else 'NO':>10s}"
                        for m in mx.METRICS)
        print(f"{label:<34s}{cells}")

    print("\nClados en comun entre pares de metricas (de 9)")
    print(f"{'':<10s}" + "".join(f"{m:>10s}" for m in mx.METRICS))
    for x in mx.METRICS:
        print(f"{x:<10s}" + "".join(f"{len(clades[x] & clades[y]):>10d}" for y in mx.METRICS))

    base = clades["gident"]
    for m in mx.METRICS:
        if clades[m] == base:
            continue
        print(f"\ndiferencias de {m} respecto de gident")
        for c in sorted(clades[m] - base, key=len):
            print(f"  solo {m}:   {{{', '.join(sorted(c))}}}")
        for c in sorted(base - clades[m], key=len):
            print(f"  solo gident: {{{', '.join(sorted(c))}}}")

    # pident = gident and bitscore = evalue, so two panels cover the four.
    plots.tree_comparison(EXPECTED_NEWICK, {
        'Identity (pident, gident)': newicks['pident'],
        'Bit score, E-value': newicks['bitscore'],
    })
