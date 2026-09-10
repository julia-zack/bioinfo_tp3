"""Ejercicio 2: de la matriz de distancias al arbol guia."""

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import matrices as mx
import plots
import tree as guide
from blast import run_blast
from exercises.ex1 import sequences

def a():
    """Construye el arbol guia de cada metrica y lo dibuja."""
    _, ids = sequences()
    raw = mx.build(run_blast(), ids)

    trees = {}
    for metric in mx.METRICS:
        D = mx.to_distance(mx.symmetrize(raw[metric], metric), metric)
        newick = guide.build(D, ids)
        trees[metric] = newick
        print(f"\n{metric}\n{newick}")

        plots.guide_tree(newick, D, ids,
                         f"Arbol guia -- {mx.LABELS[metric]}",
                         f"ex2_tree_{metric}.png")
    return trees
