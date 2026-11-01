"""From BLAST's list of pairs to the matrices of the consigna: one per score,
their symmetry, and their conversion into distances the guide tree can use.

The three scores are not three views of one number. Identity is a property of
the alignment, bit score is a property of the alignment *and* the substitution
matrix, and the E-value is a property of the alignment, the substitution
matrix, the query length and the size of the database it was searched against.
Only the first two are candidates for a distance; the third is kept because
comparing it is the point of the exercise.
"""

import math

import numpy as np

from sequences import gident

# The metrics carried all the way through: BLAST's identity, the coverage
# corrected identity of sequences.gident, the bit score and the E-value.
METRICS = ["pident", "gident", "bitscore", "evalue"]

LABELS = {
    "pident":   "Identidad BLAST (% sobre el HSP)",
    "gident":   "Identidad global (% sobre la mas corta)",
    "bitscore": "Bit score",
    "evalue":   "E-value",
}

# What an absent pair means, per metric. BLAST reporting no hit is not the
# same as BLAST reporting a bad one, and the difference has to be recorded
# somewhere rather than silently becoming a zero.
# The same metrics, named for figure titles.
NAMES = {
    "pident":   "BLAST identity (% over the HSP)",
    "gident":   "Global identity (% over the shorter sequence)",
    "bitscore": "Bit score",
    "evalue":   "E-value",
}

MISSING = {"pident": 0.0, "gident": 0.0, "bitscore": 0.0, "evalue": np.inf}


def build(hits, ids):
    """{metric: 10x10 array}, indexed by the order of `ids`.

    The arrays are the *raw* asymmetric ones: entry [i, j] is what BLAST
    reported using i as query and j as subject, and [j, i] is the other run.
    Nothing is averaged here on purpose, since whether those two differ is
    itself one of the questions.
    """
    index = {name: k for k, name in enumerate(ids)}
    n = len(ids)

    out = {m: np.full((n, n), MISSING[m], dtype=float) for m in METRICS}
    for row in hits:
        i, j = index[row["qseqid"]], index[row["sseqid"]]
        out["pident"][i, j] = row["pident"]
        out["gident"][i, j] = gident(row)
        out["bitscore"][i, j] = row["bitscore"]
        out["evalue"][i, j] = row["evalue"]
    return out


def missing_pairs(hits, ids):
    """The ordered pairs BLAST did not align at all."""
    seen = {(h["qseqid"], h["sseqid"]) for h in hits}
    return [(a, b) for a in ids for b in ids if a != b and (a, b) not in seen]


# ---------------------------------------------------------------------------
# Symmetry
# ---------------------------------------------------------------------------

def symmetry(M):
    """How far a raw matrix is from symmetric: (n_diff, max_abs, mean_rel).

    Measured only off the diagonal and only where both directions exist.
    `mean_rel` is relative because the three metrics live on different scales:
    a bit score off by 2 and an E-value off by 2 are not comparable failures.
    """
    n = M.shape[0]
    diffs, rels = [], []
    for i in range(n):
        for j in range(i + 1, n):
            a, b = M[i, j], M[j, i]
            if not (np.isfinite(a) and np.isfinite(b)):
                continue
            diffs.append(abs(a - b))
            scale = max(abs(a), abs(b))
            if scale > 0:
                rels.append(abs(a - b) / scale)
    if not diffs:
        return 0, float("nan"), float("nan")
    n_diff = sum(1 for d in diffs if d > 0)
    return n_diff, max(diffs), (sum(rels) / len(rels) if rels else 0.0)


def symmetrize(M, metric):
    """The average of the two directions, which is what the tree gets fed.

    Averaging, rather than picking one direction, because neither direction is
    the right one: 'the distance from human to chicken' does not have a query
    and a subject. The E-value is averaged in log space, since it is a
    probability spanning 200 orders of magnitude and (1e-100 + 1e-3)/2 is just
    1e-3 with extra steps.
    """
    S = M.copy()
    n = M.shape[0]
    for i in range(n):
        for j in range(i + 1, n):
            a, b = M[i, j], M[j, i]
            if metric == "evalue":
                avg = 10 ** ((safe_log10(a) + safe_log10(b)) / 2)
            else:
                avg = (a + b) / 2
            S[i, j] = S[j, i] = avg
    return S


LOG_FLOOR = 1e-200   # blastp prints 0.0 for a perfect self-hit


def safe_log10(e):
    return math.log10(max(float(e), LOG_FLOOR)) if np.isfinite(e) else 0.0


# ---------------------------------------------------------------------------
# Score -> distance
# ---------------------------------------------------------------------------

def to_distance(S, metric):
    """A symmetric matrix of similarities into one of distances in [0, 1].

    Each metric needs its own conversion, and getting this wrong is the
    quietest way to end up with a beautiful, backwards tree:

    identity   d = 1 - id. The textbook p-distance. It saturates: near the
               divergence of two random sequences (~95% for 20 amino acids;
               the usual 75% is the DNA value) every pair looks equally far,
               which is why real phylogenetics applies a Poisson or Kimura
               correction on top.

    bit score  not a similarity in [0, 1] at all, it grows with length, so a
               long pair of distant sequences can outscore a short pair of
               near-identical ones. Dividing by the smaller of the two
               self-scores normalises that away: the self-score is the most
               the pair could possibly have scored.

    E-value    the awkward one. -log10(E) is monotone in similarity, so it can
               be rescaled into [0, 1], but the result is not a distance in
               any useful sense: it depends on the size of the database, so
               the same two sequences change distance when a sequence is added
               to the set, and it tends to 0.0 for anything very close.
               The scaling below makes the matrix comparable to the others; the report should say it
               is a ranking, not a metric.
    """
    n = S.shape[0]
    D = np.zeros((n, n))

    if metric in ("pident", "gident"):
        D = 1.0 - S / 100.0

    elif metric == "bitscore":
        self_score = np.diag(S).copy()
        for i in range(n):
            for j in range(n):
                ref = min(self_score[i], self_score[j])
                D[i, j] = 1.0 - S[i, j] / ref if ref > 0 else 1.0

    elif metric == "evalue":
        logs = np.array([[safe_log10(S[i, j]) for j in range(n)] for i in range(n)])
        off = [logs[i, j] for i in range(n) for j in range(n) if i != j
               and np.isfinite(logs[i, j])]
        lo, hi = min(off), max(off)          # hi is the least significant pair
        span = (hi - lo) or 1.0
        D = (logs - lo) / span               # more significant -> smaller
        np.fill_diagonal(D, 0.0)

    D = np.clip(D, 0.0, 1.0)
    np.fill_diagonal(D, 0.0)
    return D


# ---------------------------------------------------------------------------
# Comparing matrices
# ---------------------------------------------------------------------------

def upper(M):
    """The off-diagonal half as a flat vector: the 45 numbers that actually
    determine the tree."""
    n = M.shape[0]
    return np.array([M[i, j] for i in range(n) for j in range(i + 1, n)])


def spearman(x, y):
    """Rank correlation, written out because scipy is not a dependency here.

    Rank rather than Pearson because the tree never looks at the values, only
    at which pair is smallest. Two matrices that disagree on every number but
    agree on the ordering build the identical tree.
    """
    rx, ry = rank(x), rank(y)
    rx = rx - rx.mean()
    ry = ry - ry.mean()
    denom = math.sqrt((rx ** 2).sum() * (ry ** 2).sum())
    return float((rx * ry).sum() / denom) if denom else float("nan")


def rank(v):
    """Ranks with ties averaged."""
    order = np.argsort(v, kind="mergesort")
    ranks = np.empty(len(v), dtype=float)
    ranks[order] = np.arange(len(v), dtype=float)
    # average the ranks of tied values
    _, inverse, counts = np.unique(v, return_inverse=True, return_counts=True)
    sums = np.zeros(len(counts))
    np.add.at(sums, inverse, ranks)
    return (sums / counts)[inverse]


def nearest_neighbours(D, ids):
    """For each sequence, its closest partner. This is the part of the matrix
    the first step of the guide tree consumes, so two matrices that produce
    the same table produce the same first fusion."""
    out = {}
    for i, name in enumerate(ids):
        order = [j for j in np.argsort(D[i]) if j != i]
        out[name] = ids[order[0]]
    return out


# ---------------------------------------------------------------------------
# Printing
# ---------------------------------------------------------------------------

def show(M, ids, title, fmt="{:8.3f}", half=False):
    """The matrix, as the consigna draws it: names down the side, and
    optionally only the relevant half below the diagonal."""
    width = max(len(s) for s in ids)
    print(f"\n{title}")
    print(" " * (width + 1) + "".join(f"{s[:8]:>9s}" for s in ids))
    for i, name in enumerate(ids):
        cells = []
        for j in range(len(ids)):
            if half and j >= i:
                cells.append(" " * 9)
            elif np.isinf(M[i, j]):
                cells.append(f"{'inf':>9s}")
            else:
                cells.append(fmt.format(M[i, j]).rjust(9))
        print(f"{name:<{width}s} " + "".join(cells))


# ---------------------------------------------------------------------------
# From the distance matrix to the guide tree
#
# From the distance matrix to the guide tree: the procedure of the assignment.
#
# Repeatedly:
#   i)  find the two closest elements and join them in a node;
#   ii) recompute the distance from the "other" elements to the new node as the
#       average of their distances to its two children, and look for the closest
#       pair again.
#
# The two cases the assignment distinguishes -- option a) the two closest are
# loose elements, option b) they are an element and an existing node -- need no
# separate branches in the code: an existing node IS an element of the new
# matrix, so "find the closest pair" covers both.
# ---------------------------------------------------------------------------



class Cluster:
    """An element of the matrix: first a sequence, later a node.

    `height` is the node's height, half the distance at which its two
    children were joined. Leaves are at height 0.
    """

    def __init__(self, newick, height=0.0, leaves=None, children=None):
        self.newick = newick
        self.height = height
        self.leaves = leaves or {newick}
        # A leaf is the name of the sequence; a node is the tuple of its two
        # children. It is the same tree as the Newick one, but walkable:
        # exercise 3 needs it to know in which order to merge.
        self.nested = children if children is not None else newick


def newick(D, ids):
    """The guide tree in Newick format, from a symmetric matrix."""
    return build_tree(D, ids).newick + ";"


def build_tree(D, ids):
    """The root node of the guide tree: `.newick` to draw it, `.nested` to
    walk it (leaves = names, nodes = (child_a, child_b))."""
    n = len(ids)
    active = list(range(n))
    clusters = {i: Cluster(ids[i]) for i in range(n)}
    dist = {frozenset((i, j)): float(D[i][j])
            for i in range(n) for j in range(i + 1, n)}
    next_id = n

    while len(active) > 1:
        # i) the closest pair. Ties are broken by name so that two runs on the
        # same matrix give the same tree: with identical sequences (human and
        # chimp) the distance is 0 and the tie is real.
        pair = min(
            (frozenset((a, b)) for k, a in enumerate(active) for b in active[k + 1:]),
            key=lambda p: (dist[p], sorted(n for k in p for n in clusters[k].leaves)),
        )
        a, b = sorted(pair, key=lambda x: sorted(clusters[x].leaves))
        height = dist[pair] / 2

        # The branch to each child is what is missing from the child's height
        # up to the parent's. Averaging never gives a distance below the pair
        # just merged, so heights only grow; the max(..., 0) is a safeguard.
        la = max(height - clusters[a].height, 0.0)
        lb = max(height - clusters[b].height, 0.0)

        clusters[next_id] = Cluster(
            newick=f"({clusters[a].newick}:{la:.4f},{clusters[b].newick}:{lb:.4f})",
            height=height,
            leaves=clusters[a].leaves | clusters[b].leaves,
            children=(clusters[a].nested, clusters[b].nested),
        )

        # ii) recompute the distance from the new node to the rest.
        active.remove(a)
        active.remove(b)
        for other in active:
            dist[frozenset((next_id, other))] = (
                dist[frozenset((a, other))] + dist[frozenset((b, other))]) / 2
        active.append(next_id)
        next_id += 1

    return clusters[active[0]]


def clades(nested):
    """The set of clades (each one a frozenset of leaves) of a nested tree.

    This is what two trees are compared by: the order of the two children of
    a node carries no information.
    """
    found = set()

    def leaves(node):
        if isinstance(node, str):
            return frozenset([node])
        below = leaves(node[0]) | leaves(node[1])
        found.add(below)
        return below

    leaves(nested)
    return found
