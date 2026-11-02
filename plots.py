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
MSA_COLUMNS_PER_ROW = 50

DISTANCE_CMAP = LinearSegmentedColormap.from_list('distance', ['#ffffff', BLUE])

HERE = os.path.dirname(os.path.abspath(__file__))
FIGURES = os.path.join(HERE, 'figures')


def _tidy(ax):
    ax.yaxis.grid(True, color=GRID, linewidth=0.8)
    ax.set_axisbelow(True)
    for side in ('top', 'right'):
        ax.spines[side].set_visible(False)


def _save(fig, filename):
    os.makedirs(FIGURES, exist_ok=True)
    fig.savefig(os.path.join(FIGURES, filename), dpi=150, bbox_inches='tight')


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
        ax.set_xlabel(f'distance ({reference})')
        ax.set_title(metric, fontsize=10)
        ax.set_xlim(-0.02, 1.02)
        ax.set_ylim(-0.02, 1.02)
        _tidy(ax)

    axes[0].set_ylabel('distance (panel metric)')
    fig.suptitle('Agreement between metrics, one point per sequence pair')
    fig.tight_layout()
    _save(fig, filename)
    plt.show()


def guide_tree(newick, D, ids, title, filename='ex2_guide_tree.png'):
    """The guide tree and, next to it, the distance matrix it came from.

    Phylo draws the branches to scale, so the horizontal position of a node is
    its accumulated distance from the root.

    The rows of the matrix are reordered to follow the order of the tree's
    leaves, which is what makes putting them side by side mean something: each
    row sits at the height of its leaf, and the light blocks along the
    diagonal match the clades of the tree.
    """
    from io import StringIO
    from Bio import Phylo

    t = Phylo.read(StringIO(newick), 'newick')
    t.ladderize()   # sorts clades by size: same tree, easier to read

    # Top to bottom in the drawing, which is the order Phylo walks the leaves.
    order = [leaf.name for leaf in t.get_terminals()]
    index = [ids.index(name) for name in order]
    M = np.asarray(D)[np.ix_(index, index)]

    fig, (ax_tree, ax_mat) = plt.subplots(
        1, 2, figsize=(15, 6), gridspec_kw={'width_ratios': [1.1, 1]})

    Phylo.draw(t, axes=ax_tree, do_show=False,
               # internal nodes have no name and Phylo would write "Clade"
               # on each one
               label_func=lambda c: c.name if c.is_terminal() else "")
    ax_tree.set_xlabel('distance')
    ax_tree.set_ylabel('')
    ax_tree.set_yticks([])
    for side in ('top', 'right', 'left'):
        ax_tree.spines[side].set_visible(False)
    ax_tree.xaxis.grid(True, color=GRID, linewidth=0.8)
    ax_tree.set_axisbelow(True)

    top = float(np.max(M))
    im = ax_mat.imshow(M, cmap=DISTANCE_CMAP, vmin=0, vmax=top)

    # The value in each cell, without the leading zero (".007", not "0.007"):
    # every distance is in [0, 1], so that zero tells nothing and only takes
    # the width the third decimal needs.
    for r in range(len(order)):
        for c in range(len(order)):
            if r == c:
                continue          # la diagonal es 0 por construccion
            value = M[r, c]
            ax_mat.text(c, r, f"{value:.3f}".lstrip("0"),
                        ha='center', va='center', fontsize=6,
                        # white on the dark cells, so it stays readable
                        color='white' if value > 0.6 * top else '#333333')

    ax_mat.set_xticks(range(len(order)))
    ax_mat.set_yticks(range(len(order)))
    ax_mat.set_xticklabels([s[:12] for s in order], rotation=90, fontsize=7)
    ax_mat.set_yticklabels([s[:12] for s in order], fontsize=7)
    for side in ax_mat.spines.values():
        side.set_visible(False)
    fig.colorbar(im, ax=ax_mat, fraction=0.046, pad=0.04, label='distance')

    fig.suptitle(title)
    fig.tight_layout()
    _save(fig, filename)
    plt.show()


def tree_comparison(expected, obtained, filename='ex2_expected_vs_obtained.png'):
    """The expected tree next to the obtained ones, topology only.

    `expected` is a Newick string and `obtained` maps a panel title to one.
    Branch lengths are dropped (the expected tree has none), so every branch
    has the same length and only who joins with whom is drawn. In the obtained
    trees, the branch of every group that is not in the expected tree is
    orange.
    """
    from io import StringIO
    from Bio import Phylo

    def read(newick):
        t = Phylo.read(StringIO(newick), 'newick')
        t.ladderize()
        for c in t.find_clades():
            c.branch_length = None
        return t

    def groups(t):
        return {frozenset(leaf.name for leaf in c.get_terminals())
                for c in t.find_clades()}

    reference = read(expected)
    good = groups(reference)

    panels = [('Expected', reference)] + [(title, read(n)) for title, n in obtained.items()]
    fig, axes = plt.subplots(1, len(panels), figsize=(5.5 * len(panels), 5))

    for ax, (title, t) in zip(axes, panels):
        # Every clade gets its colour explicitly: Phylo passes a clade's
        # colour down to its children, so colouring only the wrong groups
        # would paint everything below them too.
        for c in t.find_clades():
            leaves = frozenset(leaf.name for leaf in c.get_terminals())
            c.color = 'black' if leaves in good else ORANGE
        Phylo.draw(t, axes=ax, do_show=False,
                   label_func=lambda c: c.name if c.is_terminal() else "")
        ax.set_title(title, fontsize=11)
        ax.set_xlabel('')
        ax.set_ylabel('')
        ax.set_xticks([])
        ax.set_yticks([])
        for side in ax.spines.values():
            side.set_visible(False)

    fig.suptitle('Expected vs obtained guide trees (orange: group not in the expected tree)')
    fig.tight_layout()
    _save(fig, filename)
    plt.show()


def msa(alignment, title, filename='ex3_msa.png', per_row=MSA_COLUMNS_PER_ROW,
        symbols=None):
    """The MSA in blocks of `per_row` columns, one letter per cell.

    The background of each cell is the fraction of the column that shares its
    residue (white = nobody else has it, blue = the whole column): a column
    that is blue from top to bottom is a conserved position, readable without
    reading the letters. Gaps are grey. `symbols`, one character per column
    (the Clustal '*', ':', '.'), goes in an extra line under each block.
    """
    names = list(alignment)
    n_cols = len(alignment[names[0]])
    n_blocks = -(-n_cols // per_row)
    n_lines = len(names) + (symbols is not None)

    fig, axes = plt.subplots(n_blocks, 1, figsize=(per_row * 0.24 + 2.2, n_blocks * (n_lines * 0.27 + 0.55)))
    axes = np.atleast_1d(axes)

    for block, ax in enumerate(axes):
        start = block * per_row
        stop = min(start + per_row, n_cols)
        shade = np.ones((len(names), stop - start))
        for c in range(start, stop):
            column = [alignment[n][c] for n in names]
            for r, ch in enumerate(column):
                if ch == '-':
                    shade[r, c - start] = np.nan
                else:
                    shade[r, c - start] = column.count(ch) / len(names)
        cmap = DISTANCE_CMAP.copy()
        cmap.set_bad(GRID)
        ax.imshow(np.ma.masked_invalid(shade), cmap=cmap, vmin=0, vmax=1, aspect='auto')
        for r, name in enumerate(names):
            for c in range(start, stop):
                ch = alignment[name][c]
                ax.text(c - start, r, ch, ha='center', va='center', fontsize=7.5,
                        family='monospace',
                        color='#999999' if ch == '-' else
                        ('white' if shade[r, c - start] > 0.75 else '#333333'))
        if symbols is not None:
            for c in range(start, stop):
                ax.text(c - start, len(names), symbols[c], ha='center', va='center',
                        fontsize=9, family='monospace', fontweight='bold', color='#333333')
            ax.set_ylim(len(names) + 0.5, -0.5)
        ax.set_yticks(range(len(names)))
        ax.set_yticklabels(names, fontsize=7.5)
        ax.set_xticks(range(0, stop - start, 10))
        ax.set_xticklabels([str(start + k + 1) for k in range(0, stop - start, 10)],
                           fontsize=7)
        ax.tick_params(length=0)
        for side in ax.spines.values():
            side.set_visible(False)

    fig.suptitle(title)
    fig.tight_layout()
    _save(fig, filename)
    plt.show()


def _upper(M):
    n = M.shape[0]
    return np.array([M[i, j] for i in range(n) for j in range(i + 1, n)])
