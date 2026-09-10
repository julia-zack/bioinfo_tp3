"""Plot helpers reused across exercises.

Same conventions as TP1: a recessive grid, no top/right spines, and the same
two accent colours, so the figures of both reports look like one series.
"""

import os

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

BLUE = '#4a6fa5'
ORANGE = '#d1852f'
GREY = '#b0b0b0'
GRID = '#e5e5e5'

# White (identical) to blue (unrelated). Sequential rather than a rainbow: the
# quantity has an order, so the colour has to have one too, and viridis-style
# hue rotation makes near-equal cells look categorically different.
DISTANCE_CMAP = LinearSegmentedColormap.from_list('distance', ['#ffffff', BLUE])

HERE = os.path.dirname(os.path.abspath(__file__))
FIGURES = os.path.join(HERE, 'exercises')


def _tidy(ax):
    ax.yaxis.grid(True, color=GRID, linewidth=0.8)
    ax.set_axisbelow(True)
    for side in ('top', 'right'):
        ax.spines[side].set_visible(False)


def _save(fig, filename):
    fig.savefig(os.path.join(FIGURES, filename), dpi=150, bbox_inches='tight')


def heatmaps(distances, ids, filename='ex1_distance_matrices.png'):
    """One panel per metric, on a shared 0-1 scale.

    Shared scale on purpose: the point of the figure is that the same 10
    sequences look differently spread depending on the score, and per-panel
    normalisation would hide exactly that by stretching each one to fill the
    colour range.
    """
    metrics = list(distances)
    fig, axes = plt.subplots(1, len(metrics), figsize=(4.2 * len(metrics), 4.6))
    axes = np.atleast_1d(axes)

    for k, (ax, metric) in enumerate(zip(axes, metrics)):
        D = distances[metric]
        im = ax.imshow(D, cmap=DISTANCE_CMAP, vmin=0, vmax=1)
        ax.set_title(metric, fontsize=10)
        ax.set_xticks(range(len(ids)))
        ax.set_yticks(range(len(ids)))
        ax.set_xticklabels([s[:12] for s in ids], rotation=90, fontsize=7)
        # Only the first panel is labelled: the rows are the same 10 sequences
        # in the same order everywhere, and repeating them four times pushes
        # each label over the panel to its left.
        ax.set_yticklabels([s[:12] for s in ids] if k == 0 else [], fontsize=7)
        for side in ax.spines.values():
            side.set_visible(False)

    fig.colorbar(im, ax=axes, fraction=0.02, pad=0.02, label='distancia')
    fig.suptitle('Matrices de distancia segun el puntaje de BLAST')
    _save(fig, filename)
    plt.show()


def scatter(distances, ids, reference='gident',
            filename='ex1_metric_agreement.png'):
    """Every other metric against the reference, one point per pair.

    The diagonal line is where the two metrics agree. What matters is not the
    scatter around it but whether the *order* of the points is preserved:
    points that swap position are pairs whose fusion order in the guide tree
    changes.
    """
    others = [m for m in distances if m != reference]
    ref = _upper(distances[reference])

    fig, axes = plt.subplots(1, len(others), figsize=(4.2 * len(others), 4),
                             sharey=True)
    axes = np.atleast_1d(axes)

    for ax, metric in zip(axes, others):
        y = _upper(distances[metric])
        ax.plot([0, 1], [0, 1], color=GREY, linestyle='--', linewidth=1, zorder=0)
        ax.scatter(ref, y, color=BLUE, alpha=0.75,
                   edgecolors='white', linewidths=0.5, zorder=3)
        ax.set_xlabel(f'distancia ({reference})')
        ax.set_title(metric, fontsize=10)
        ax.set_xlim(-0.02, 1.02)
        ax.set_ylim(-0.02, 1.02)
        _tidy(ax)

    axes[0].set_ylabel('distancia (metrica del panel)')
    fig.suptitle('Acuerdo entre metricas, un punto por par de secuencias')
    fig.tight_layout()
    _save(fig, filename)
    plt.show()


def guide_tree(newick, D, ids, title, filename='ex2_guide_tree.png'):
    """El arbol guia y, al lado, la matriz de distancias de la que salio.

    Phylo dibuja las ramas a escala de sus largos, asi que la posicion
    horizontal de cada nodo es su distancia acumulada desde la raiz.

    Las filas de la matriz se reordenan siguiendo el orden de las hojas del
    arbol, que es lo que hace que ponerlas al lado signifique algo: cada fila
    queda a la altura de su hoja, y los bloques claros sobre la diagonal se
    corresponden con los clados del arbol.
    """
    from io import StringIO
    from Bio import Phylo

    t = Phylo.read(StringIO(newick), 'newick')
    t.ladderize()   # ordena los clados por tamano: mismo arbol, mas legible

    # De arriba hacia abajo en el dibujo, que es el orden en que Phylo
    # recorre las hojas.
    order = [leaf.name for leaf in t.get_terminals()]
    index = [ids.index(name) for name in order]
    M = np.asarray(D)[np.ix_(index, index)]

    fig, (ax_tree, ax_mat) = plt.subplots(
        1, 2, figsize=(15, 6), gridspec_kw={'width_ratios': [1.1, 1]})

    Phylo.draw(t, axes=ax_tree, do_show=False,
               # los nodos internos no tienen nombre y Phylo escribiria
               # "Clade" en cada uno
               label_func=lambda c: c.name if c.is_terminal() else "")
    ax_tree.set_xlabel('distancia')
    ax_tree.set_ylabel('')
    ax_tree.set_yticks([])
    for side in ('top', 'right', 'left'):
        ax_tree.spines[side].set_visible(False)
    ax_tree.xaxis.grid(True, color=GRID, linewidth=0.8)
    ax_tree.set_axisbelow(True)

    top = float(np.max(M))
    im = ax_mat.imshow(M, cmap=DISTANCE_CMAP, vmin=0, vmax=top)

    # El valor en cada celda. Sin el cero inicial (".007" en vez de "0.007")
    # porque todas las distancias estan en [0, 1] y ese cero no distingue
    # nada, solo ocupa el ancho que necesita el tercer decimal.
    for r in range(len(order)):
        for c in range(len(order)):
            if r == c:
                continue          # la diagonal es 0 por construccion
            value = M[r, c]
            ax_mat.text(c, r, f"{value:.3f}".lstrip("0"),
                        ha='center', va='center', fontsize=6,
                        # blanco sobre las celdas oscuras, para que se lea
                        color='white' if value > 0.6 * top else '#333333')

    ax_mat.set_xticks(range(len(order)))
    ax_mat.set_yticks(range(len(order)))
    ax_mat.set_xticklabels([s[:12] for s in order], rotation=90, fontsize=7)
    ax_mat.set_yticklabels([s[:12] for s in order], fontsize=7)
    for side in ax_mat.spines.values():
        side.set_visible(False)
    fig.colorbar(im, ax=ax_mat, fraction=0.046, pad=0.04, label='distancia')

    fig.suptitle(title)
    fig.tight_layout()
    _save(fig, filename)
    plt.show()


def _upper(M):
    n = M.shape[0]
    return np.array([M[i, j] for i in range(n) for j in range(i + 1, n)])
