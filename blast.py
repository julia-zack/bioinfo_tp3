"""The all-vs-all comparison: runs blastp of the family against itself and
parses the result into one row per ordered pair.

Needs BLAST+ installed and on the PATH (`brew install blast` on macOS,
`apt install ncbi-blast+` on Debian/Ubuntu). Check with `blastp -version`.
"""

import os
import json
import subprocess

from ncbi import DATA_DIR, FASTA_CACHE

HITS_CACHE = os.path.join(DATA_DIR, 'blast_hits.json')

# The columns asked of `-outfmt 6`, in order. qlen/slen and nident are not in
# the default set and are the ones that make the identity question answerable:
# see gident() below.
FIELDS = ["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
          "qstart", "qend", "sstart", "send", "evalue", "bitscore",
          "nident", "qlen", "slen"]

INT_FIELDS = {"length", "mismatch", "gapopen", "qstart", "qend",
              "sstart", "send", "nident", "qlen", "slen"}
FLOAT_FIELDS = {"pident", "evalue", "bitscore"}


def run_blast(fasta=FASTA_CACHE, force=False):
    """Every ordered pair (i, j) that BLAST can align, as a list of dicts.

    One `-subject` run rather than makeblastdb + one query at a time: with a
    query file of 10 and a subject file of the same 10, BLAST does all 100
    ordered comparisons in a single call, self-hits included. Keeping the
    self-hits is deliberate, they are the diagonal the bit score has to be
    normalised against (a bit score of 200 means nothing until you know the
    sequence scores 210 against itself).

    Two flags that change the answer:

    -evalue 10000  the default cutoff is 10, which quietly drops the most
                   distant pairs. A missing pair is not "distance infinity",
                   it is a hole in the matrix, and the clustering of exercise
                   2 has to invent a number to fill it. Raising the cutoff
                   moves that invention from the code to BLAST's statistics.

    -max_hsps 1    without it a pair with two local alignments contributes two
                   rows and the parser has to decide which one is "the"
                   comparison. One HSP per pair, the best one, is the same
                   convention the pairwise identity of a phylogeny assumes.
    """
    if os.path.exists(HITS_CACHE) and not force:
        print(f"Using cached BLAST hits from {HITS_CACHE}")
        with open(HITS_CACHE) as f:
            return json.load(f)

    cmd = [
        "blastp",
        "-query", fasta,
        "-subject", fasta,
        "-outfmt", "6 " + " ".join(FIELDS),
        "-evalue", "10000",
        "-max_hsps", "1",
    ]
    print("$ " + " ".join(cmd))
    out = subprocess.run(cmd, capture_output=True, text=True, check=True).stdout

    hits = [parse_row(line) for line in out.strip().split("\n") if line]
    os.makedirs(DATA_DIR, exist_ok=True)
    with open(HITS_CACHE, "w") as f:
        json.dump(hits, f, indent=2)
    print(f"{len(hits)} ordered pairs aligned")
    return hits


def parse_row(line):
    row = dict(zip(FIELDS, line.split("\t")))
    for k in INT_FIELDS:
        row[k] = int(row[k])
    for k in FLOAT_FIELDS:
        row[k] = float(row[k])
    return row


def gident(row):
    """Identity over the shorter sequence, not over the aligned region.

    BLAST's own `pident` is nident/length, where length is the length of the
    HSP. It answers "of the part that aligned, how much matched?", which is
    the wrong question for a distance: two sequences sharing one conserved
    30-residue motif and nothing else report pident 95 and land next to each
    other in the tree. Dividing by the shorter sequence instead folds coverage
    back in, so a partial alignment cannot buy a short distance.

    Both numbers are built into matrices, precisely so the report can show
    what the difference does to the tree.
    """
    return 100.0 * row["nident"] / min(row["qlen"], row["slen"])
