"""From the guide tree to the multiple alignment (MSA), using the pairwise
alignments (APs) BLAST produced in exercise 1.

The MSA is built by pasting pairwise alignments together in the order the
tree gives, and nothing is re-aligned. The three rules of the assignment:

  i)   two sequences that share a node enter the MSA as in their AP;
  ii)  a sequence joining a node that is already aligned keeps its AP with
       *one* of the node's sequences (the one with the best score);
  iii) two nodes are joined by the best AP between any pair (one sequence from
       each side), without touching the MSA inside each node.

They are one single case. A block is a set of rows that are already aligned (a
lone sequence is a block of one row), and `merge` joins two blocks A and B by
the AP of some a in A with some b in B. Rule i) is A and B with one row each,
rule ii) is B with one row.

One thing the assignment does not settle: BLAST is a *local* aligner, so its AP
can leave residues at the ends unaligned, and the MSA has to place all of them.
See `full_length`.
"""

from Bio.Align import substitution_matrices

GAP = "-"
SP_GAP = -4          # cost of a residue-gap pair in the sum-of-pairs score
LINE_WIDTH = 60      # columns per block when the MSA is printed


# ---------------------------------------------------------------------------
# From BLAST's row to the full-length AP
# ---------------------------------------------------------------------------

def full_length(hit, seqs):
    """BLAST's HSP, extended to cover both sequences entirely.

    BLAST reports only the local stretch, e.g. with `qstart=5` the first 4
    residues of the query are not in the alignment. An MSA contains every
    residue of every sequence, so the ends BLAST left out have to go somewhere.
    What we do is what a global alignment with free end gaps would do: the ends
    of the two sequences are laid against each other, flush with the HSP, and
    what is left of the longer one goes against gaps.

        zebrafish   MVEW|TDAERT...
        cow         --ML|TAEEKA...

    Laying the ends against each other is our decision, not BLAST's: they are
    not part of the HSP and no score backs them. They are few residues (at most
    7 in total in this data).
    """
    q, s = hit["qseqid"], hit["sseqid"]
    qfull, sfull = seqs[q], seqs[s]

    q_left, s_left = qfull[:hit["qstart"] - 1], sfull[:hit["sstart"] - 1]
    q_right, s_right = qfull[hit["qend"]:], sfull[hit["send"]:]

    left = max(len(q_left), len(s_left))
    right = max(len(q_right), len(s_right))
    row_q = q_left.rjust(left, GAP) + hit["qseq"] + q_right.ljust(right, GAP)
    row_s = s_left.rjust(left, GAP) + hit["sseq"] + s_right.ljust(right, GAP)

    assert len(row_q) == len(row_s)
    assert row_q.replace(GAP, "") == qfull, f"{q} vs {s}: the AP lost residues of the query"
    assert row_s.replace(GAP, "") == sfull, f"{q} vs {s}: the AP lost residues of the subject"
    return row_q, row_s


def pairwise_table(hits, seqs):
    """{(a, b): (row_a, row_b, bitscore)} for every ordered pair.

    BLAST leaves two alignments per pair of sequences, a-against-b and
    b-against-a, and they do not always have the same score (exercise 1,
    symmetry: the bit score differs in 9 of the 45 pairs). The MSA needs *one*
    AP per pair: we keep the one with the higher bit score and, on a tie, the
    one whose query comes first alphabetically. Both directions point to that
    same AP, so `table[a, b]` and `table[b, a]` are the same alignment, each
    with the row of its own first sequence first.
    """
    by_pair = {(h["qseqid"], h["sseqid"]): h for h in hits
               if h["qseqid"] != h["sseqid"]}

    table = {}
    for (a, b) in by_pair:
        # each pair once, even if BLAST left only one direction
        if a > b and (b, a) in by_pair:
            continue
        candidates = [by_pair[a, b]] + ([by_pair[b, a]] if (b, a) in by_pair else [])
        best = min(candidates, key=lambda h: (-h["bitscore"], h["qseqid"]))
        row_q, row_s = full_length(best, seqs)
        score = best["bitscore"]
        if best["qseqid"] == a:
            table[a, b] = (row_q, row_s, score)
            table[b, a] = (row_s, row_q, score)
        else:
            table[a, b] = (row_s, row_q, score)
            table[b, a] = (row_q, row_s, score)
    return table


# ---------------------------------------------------------------------------
# Joining two blocks
# ---------------------------------------------------------------------------

def merge(A, B, table):
    """Join two blocks of already aligned rows by the best AP between them.

    A and B are {name: row}; the rows of a block all have the same length.
    Returns (new_block, (a, b, score)), with (a, b) the chosen pair.

    The anchor pair (a in A, b in B) is the bridge: `a` shows up in two places,
    its row in A and the AP, and `b` in its row in B and the AP. The residues
    of `a` are the same sequence in both, so walking the row of `a` in A and
    the row of `a` in the AP side by side tells, column by column, which column
    of A matches which column of the AP (and the same for `b` and B). Three
    things are walked at once, one pointer on each, and every step emits one
    column of the MSA:

      1. A has a gap in `a`: a column of A the AP knows nothing about. It goes
         out as it is, with gaps in every row of B.
      2. B has a gap in `b`: the same for B.
      3. The AP has a gap in `a` (and a residue of `b`): the AP inserted a
         residue of B that is not aligned with anything in A. B's column goes
         out, with gaps in A.
      4. The AP has a gap in `b` (and a residue of `a`): the same with A.
      5. The AP has a residue in both: the column of A and the column of B go
         out together, in the same MSA column.

    The gaps already inside A and B are never touched: this is the "once a gap,
    always a gap" of progressive aligners.
    """
    inter_block_pairs = [(x, y) for x in A for y in B]
    a, b = min(inter_block_pairs, key=lambda p: (-table[p][2], p))
    ap_a, ap_b, score = table[a, b]

    len_a, len_b, len_p = len(A[a]), len(B[b]), len(ap_a)
    out_a = {name: [] for name in A}
    out_b = {name: [] for name in B}
    i = j = k = 0       # pointers in A, in B and in the AP

    def emit_a(col):
        for name in A:
            out_a[name].append(A[name][col])

    def emit_b(col):
        for name in B:
            out_b[name].append(B[name][col])

    def gaps_a():
        for name in A:
            out_a[name].append(GAP)

    def gaps_b():
        for name in B:
            out_b[name].append(GAP)

    while i < len_a or j < len_b or k < len_p:
        if i < len_a and A[a][i] == GAP:                     # 1
            emit_a(i); gaps_b(); i += 1
        elif j < len_b and B[b][j] == GAP:                   # 2
            gaps_a(); emit_b(j); j += 1
        elif k < len_p and ap_a[k] == GAP:                   # 3
            gaps_a(); emit_b(j); j += 1; k += 1
        elif k < len_p and ap_b[k] == GAP:                   # 4
            emit_a(i); gaps_b(); i += 1; k += 1
        else:                                                # 5
            emit_a(i); emit_b(j); i += 1; j += 1; k += 1

    merged = {name: "".join(col) for name, col in {**out_a, **out_b}.items()}

    # The guarantees of the merge, checked at every step.
    assert _project(merged, A) == A, "the merge altered block A"
    assert _project(merged, B) == B, "the merge altered block B"
    assert _project(merged, {a: ap_a, b: ap_b}) == {a: ap_a, b: ap_b}, \
        "the merge does not respect the chosen AP"
    return merged, (a, b, score)


def _project(msa, block):
    """The rows of `block` inside `msa`, without the columns where every row of
    `block` has a gap. If `msa` respects `block`, the result is `block`."""
    names = list(block)
    keep = [c for c in range(len(msa[names[0]]))
            if any(msa[n][c] != GAP for n in names)]
    return {n: "".join(msa[n][c] for c in keep) for n in names}


# ---------------------------------------------------------------------------
# Walking the tree
# ---------------------------------------------------------------------------

def build(nested, seqs, table):
    """The MSA of the whole tree and the log of the merges, in order.

    `nested` is the tree of `matrices.build_tree(...).nested`: a leaf is a name, a
    node is the tuple of its two children. It is solved by recursion, children
    first, which is exactly "a node can be merged once both its children are
    aligned".
    """
    log = []

    def visit(node):
        if isinstance(node, str):
            return {node: seqs[node]}
        left, right = visit(node[0]), visit(node[1])
        merged, (a, b, score) = merge(left, right, table)
        log.append({"left": list(left), "right": list(right),
                    "anchor": (a, b), "score": score})
        return merged

    msa = visit(nested)
    check(msa, seqs)
    return msa, log


def check(msa, seqs):
    """A well-formed MSA: all rows the same length, no all-gap column, and
    removing the gaps of each row gives back its sequence."""
    lengths = {len(r) for r in msa.values()}
    assert len(lengths) == 1, f"rows of different lengths: {lengths}"
    n = lengths.pop()
    assert set(msa) == set(seqs), "missing or extra sequences"
    for name, row in msa.items():
        assert row.replace(GAP, "") == seqs[name], f"{name}: the row is not the sequence"
    for c in range(n):
        assert any(r[c] != GAP for r in msa.values()), f"column {c} is all gaps"


# ---------------------------------------------------------------------------
# Showing it and measuring it
# ---------------------------------------------------------------------------

# Clustal's similarity groups: ':' if every residue of the column falls in one
# "strong" group, '.' if it falls in one "weak" group. They are the ones the
# Clustal Omega output uses, so the two outputs read the same way.
STRONG = ["STA", "NEQK", "NHQK", "NDEQ", "QHRK", "MILV", "MILF", "HY", "FYW"]
WEAK = ["CSA", "ATV", "SAG", "STNK", "STPA", "SGND", "SNDEQK", "NDEQHK",
        "NEQHRK", "FVLIM", "HFY"]


def conservation(msa):
    """One symbol per column: '*' identical, ':' very similar, '.' somewhat
    similar, ' ' the rest (and any column with a gap)."""
    rows = list(msa.values())
    out = []
    for col in zip(*rows):
        letters = set(col)
        if GAP in letters:
            out.append(" ")
        elif len(letters) == 1:
            out.append("*")
        elif any(letters <= set(g) for g in STRONG):
            out.append(":")
        elif any(letters <= set(g) for g in WEAK):
            out.append(".")
        else:
            out.append(" ")
    return "".join(out)


def show(msa, width=LINE_WIDTH):
    """The MSA in blocks of `width` columns, with the conservation line."""
    names = list(msa)
    label = max(len(n) for n in names)
    cons = conservation(msa)
    n = len(next(iter(msa.values())))
    for start in range(0, n, width):
        for name in names:
            print(f"{name:<{label}s}  {msa[name][start:start + width]}")
        print(f"{'':<{label}s}  {cons[start:start + width]}\n")


def write_fasta(msa, path):
    with open(path, "w") as f:
        for name, row in msa.items():
            f.write(f">{name}\n{row}\n")


def read_fasta(path):
    """An aligned FASTA as {name: row}, in upper case (MAFFT writes lower
    case)."""
    out, name = {}, None
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                name = line[1:].split()[0]
                out[name] = ""
            elif line:
                out[name] += line.upper()
    return out


def aligned_pairs(msa):
    """The set of residue pairs the MSA puts in the same column:
    {(name_a, pos_a, name_b, pos_b)}, with the names in alphabetical order.

    It is the unit for comparing two MSAs of the same sequences. Comparing
    column by column does not work: an MSA with one extra column at the start
    shifts every other column and looks completely different even if it says
    the same thing about which residue goes with which.
    """
    names = sorted(msa)
    pos = {n: -1 for n in names}
    pairs = set()
    for col in zip(*(msa[n] for n in names)):
        here = {}
        for n, ch in zip(names, col):
            if ch != GAP:
                pos[n] += 1
                here[n] = pos[n]
        present = sorted(here)
        for x in range(len(present)):
            for y in range(x + 1, len(present)):
                pairs.add((present[x], here[present[x]], present[y], here[present[y]]))
    return pairs


def _row_pairs(row_a, row_b):
    """The residue pairs (i, j) two aligned rows put in the same column."""
    i = j = -1
    out = set()
    for p, q in zip(row_a, row_b):
        i += p != GAP
        j += q != GAP
        if p != GAP and q != GAP:
            out.add((i, j))
    return out


def consistency(alignment, table):
    """How much of each pair's AP survives in the MSA.

    At every merge the MSA respects *exactly* one AP, the anchor pair's. The
    other 44 are respected only as far as the MSA, built from other APs,
    reproduces them. Returns (total fraction of the APs' residue pairs that are
    in the MSA, {pair: fraction}). The APs are the full-length ones, so they
    include the end residues `full_length` laid against each other.
    """
    kept = total = 0
    per_pair = {}
    for (a, b), (row_a, row_b, _) in table.items():
        if a > b:
            continue
        ap = _row_pairs(row_a, row_b)
        in_msa = _row_pairs(alignment[a], alignment[b])
        per_pair[a, b] = len(ap & in_msa) / len(ap) if ap else 1.0
        kept += len(ap & in_msa)
        total += len(ap)
    return (kept / total if total else 1.0), per_pair


def agreement(mine, other):
    """(recall, precision) of the aligned pairs of `mine` against `other`:
    what fraction of the pairs of `other` also appears in `mine`, and what
    fraction of those of `mine` appears in `other`."""
    a, b = aligned_pairs(mine), aligned_pairs(other)
    common = len(a & b)
    return (common / len(b) if b else 1.0), (common / len(a) if a else 1.0)


def sum_of_pairs(msa, matrix=None, gap=SP_GAP):
    """Sum-of-pairs score: BLOSUM62 for every pair of residues in a column,
    `gap` for every residue-gap pair and 0 for gap-gap. It is the objective an
    MSA tries to maximise; it puts two MSAs on one scale, it does not say which
    one is true."""
    matrix = matrix or substitution_matrices.load("BLOSUM62")
    total = 0
    for col in zip(*msa.values()):
        for x in range(len(col)):
            for y in range(x + 1, len(col)):
                p, q = col[x], col[y]
                if p == GAP and q == GAP:
                    continue
                total += gap if GAP in (p, q) else matrix[p, q]
    return total


def stats(msa):
    """Length, identical columns, gap characters and mean identity per pair
    (over the columns where both rows have a residue)."""
    names = list(msa)
    n = len(msa[names[0]])
    identical = 0
    for col in zip(*msa.values()):
        if GAP not in col and len(set(col)) == 1:
            identical += 1
    gaps = sum(r.count(GAP) for r in msa.values())

    ids = []
    for x in range(len(names)):
        for y in range(x + 1, len(names)):
            ra, rb = msa[names[x]], msa[names[y]]
            both = [(p, q) for p, q in zip(ra, rb) if p != GAP and q != GAP]
            if both:     # two rows with no column in common add no identity
                ids.append(sum(p == q for p, q in both) / len(both))
    return {"length": n, "identical": identical, "gaps": gaps,
            "mean_identity": sum(ids) / len(ids) if ids else 0.0}
