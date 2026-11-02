"""The sequences of the TP: fetching the protein family from NCBI (with a local
cache so the same records are not downloaded twice) and the all-vs-all BLAST
of the family against itself.

Part 1, NCBI:

The family is beta-globin, one entry per species, and the ten entries are
named explicitly below rather than searched for. Two reasons for that choice
over "the top hit of a query per species", which can quietly return a protein
that only mentions hemoglobin in its annotation (a haptoglobin, a
transcription factor):

  - the family is paralogous. Mouse alone has five reviewed beta-like globins
    (beta-1, beta-2, beta-H0, beta-H1, epsilon-Y2), rat and Xenopus and
    zebrafish have two each. The beta-globin cluster duplicated independently
    in rodents, in amphibians and in fish, so "the beta-globin of the mouse"
    is not a question a query can answer -- it is a choice, and the choice is
    the adult major chain, which UniProt annotates as such;
  - the entries are the dataset. A search re-run in six months can return
    something else; ten accessions cannot. This is the difference between a
    method that is described and a method that is reproducible.

The candidates were listed with this UniProt query, which is the one to cite
in the report:

    https://rest.uniprot.org/uniprotkb/search
      ?query=(gene:hbb OR protein_name:"hemoglobin subunit beta")
             AND reviewed:true AND (organism_id:9606 OR organism_id:9598 ...)
      &fields=accession,id,protein_name,length,organism_name
      &format=tsv

The species tree these ten are meant to recover, for exercise 2, is
(((((((human, chimp), macaque), ((mouse, rat), rabbit)), cow), chicken), frog),
fish).
"""

import os
import json
import subprocess

from Bio import Entrez, SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# NCBI requires a contact address on every request. Change it to yours.
Entrez.email = 'arielzingman@gmail.com'

HERE = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.join(HERE, 'data')

FASTA_CACHE = os.path.join(DATA_DIR, 'globins.fasta')
META_CACHE = os.path.join(DATA_DIR, 'globins.json')

# (species, UniProt accession, UniProt entry name, why this one and not its
# paralogues). Ordered roughly by distance to human, so the printed matrices
# read as a gradient and an obviously wrong tree is obvious at a glance.
FAMILY = [
    ("Homo sapiens",          "P68871", "HBB_HUMAN",  "single adult beta"),
    ("Pan troglodytes",       "P68873", "HBB_PANTR",  "single adult beta"),
    ("Macaca mulatta",        "P02026", "HBB_MACMU",  "single adult beta"),
    ("Mus musculus",          "P02088", "HBB1_MOUSE", "beta-major, adult; over beta-2/H0/H1/epsilon-Y2"),
    ("Rattus norvegicus",     "P02091", "HBB1_RAT",   "major-form, adult; over beta-2 minor"),
    ("Oryctolagus cuniculus", "P02057", "HBB_RABIT",  "beta-1/2, the only reviewed adult beta"),
    ("Bos taurus",            "P02070", "HBB_BOVIN",  "single adult beta"),
    ("Gallus gallus",         "P02112", "HBB_CHICK",  "adult beta; over the embryonic rho/epsilon"),
    ("Xenopus laevis",        "P02132", "HBB1_XENLA", "beta-major, adult; over beta-2 and the larval chain"),
    ("Danio rerio",           "Q90486", "HBB1_DANRE", "beta-A1, adult; over beta-A2"),
]

# A globin is ~146 residues. Anything outside this range means the accession
# is not what it is supposed to be (a haptoglobin is 347 aa, a transcription
# factor 373), and that should stop the run rather than flow silently into a
# distance matrix.
MIN_LEN, MAX_LEN = 120, 200


def slug(species):
    """'Homo sapiens' -> 'Homo_sapiens'. Used as the FASTA id, which is what
    BLAST echoes back in qseqid/sseqid, which is what labels the tree."""
    return species.replace(" ", "_")


def check(species, accession, record):
    """Fail loudly if the accession did not return a globin.

    Checked on the description as well as the length, because the two failure
    modes are different: a wrong accession gives the wrong protein at the
    wrong length, while a right accession on a precursor gives the right
    protein at a length that is too long.
    """
    text = record.description.lower()
    if not ("hemoglobin" in text or "globin" in text):
        raise ValueError(f"{species}: {accession} is not a globin -- "
                         f"{record.description[:70]}")
    if not (MIN_LEN <= len(record.seq) <= MAX_LEN):
        raise ValueError(f"{species}: {accession} is {len(record.seq)} aa, "
                         f"outside {MIN_LEN}-{MAX_LEN} -- {record.description[:70]}")


def fetch_family(force=False):
    """The 10 sequences, as SeqRecords with the species as id.

    Cached as data/globins.fasta. One efetch for all ten: they are named, so
    there is nothing to search for and no reason to pay ten round trips.

    Re-labelling the records here, instead of keeping NCBI's
    '>sp|P68871.2| RecName: Full=...' headers, is not cosmetic: that id
    travels through BLAST's output into the matrix row names and from there
    into the Newick string, and Newick gives '|' and ':' a meaning.
    """
    if os.path.exists(FASTA_CACHE) and not force:
        print(f"Using cached sequences from {FASTA_CACHE}")
        return list(SeqIO.parse(FASTA_CACHE, "fasta"))

    accessions = [acc for _, acc, _, _ in FAMILY]
    print(f"Downloading {len(accessions)} entries: {', '.join(accessions)}")
    handle = Entrez.efetch(db="protein", id=",".join(accessions),
                           rettype="fasta", retmode="text")
    downloaded = list(SeqIO.parse(handle, "fasta"))
    handle.close()

    # Matched by accession rather than by position: efetch usually preserves
    # the order of the id list, but "usually" is not a guarantee worth
    # betting the whole matrix on, and a silent off-by-one here would label
    # every sequence with the wrong species.
    by_accession = {}
    for record in downloaded:
        for part in record.id.split("|"):
            by_accession[part.split(".")[0]] = record

    records, meta = [], {}
    for species, accession, entry_name, reason in FAMILY:
        if accession not in by_accession:
            raise ValueError(f"{species}: NCBI returned nothing for {accession}")
        original = by_accession[accession]
        check(species, accession, original)

        records.append(SeqRecord(Seq(str(original.seq)),
                                 id=slug(species), description=original.description))
        meta[slug(species)] = {
            "species": species,
            "accession": accession,
            "entry_name": entry_name,
            "chosen_because": reason,
            "description": original.description,
            "length": len(original.seq),
        }
        print(f"  {species:24s} {accession:8s} {entry_name:12s} "
              f"{len(original.seq):4d} aa")

    os.makedirs(DATA_DIR, exist_ok=True)
    SeqIO.write(records, FASTA_CACHE, "fasta")
    with open(META_CACHE, "w") as f:
        json.dump(meta, f, indent=2)
    return records


def load_metadata():
    """What was downloaded and why that entry and not its paralogues: the
    provenance half of the answer to 'obtenga 10 secuencias de una misma
    familia'."""
    with open(META_CACHE) as f:
        return json.load(f)


# ---------------------------------------------------------------------------
# Part 2, BLAST
#
# The all-vs-all comparison: runs blastp of the family against itself and
# parses the result into one row per ordered pair.
#
# Needs BLAST+ installed and on the PATH (`brew install blast` on macOS,
# `apt install ncbi-blast+` on Debian/Ubuntu). Check with `blastp -version`.
# ---------------------------------------------------------------------------



HITS_CACHE = os.path.join(DATA_DIR, 'blast_hits.json')

# The columns asked of `-outfmt 6`, in order. qlen/slen and nident are not in
# the default set and are the ones that make the identity question answerable:
# see gident() below.
FIELDS = ["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
          "qstart", "qend", "sstart", "send", "evalue", "bitscore",
          "nident", "qlen", "slen", "qseq", "sseq"]

INT_FIELDS = {"length", "mismatch", "gapopen", "qstart", "qend",
              "sstart", "send", "nident", "qlen", "slen"}
FLOAT_FIELDS = {"pident", "evalue", "bitscore"}


def run_blast(fasta=FASTA_CACHE, force=False, matrix=None):
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

    `matrix` is the substitution matrix. None is blastp's default (BLOSUM62)
    and the cache every other exercise reads; any other name gets its own
    cache, so comparing matrices never overwrites the main run.
    """
    cache = HITS_CACHE if matrix is None else HITS_CACHE.replace(".json", f"_{matrix.lower()}.json")
    if os.path.exists(cache) and not force:
        print(f"Using cached BLAST hits from {cache}")
        with open(cache) as f:
            return json.load(f)

    cmd = [
        "blastp",
        "-query", fasta,
        "-subject", fasta,
        "-outfmt", "6 " + " ".join(FIELDS),
        "-evalue", "10000",
        "-max_hsps", "1",
    ]
    if matrix is not None:
        cmd += ["-matrix", matrix]
    print("$ " + " ".join(cmd))
    out = subprocess.run(cmd, capture_output=True, text=True, check=True).stdout

    hits = [parse_row(line) for line in out.strip().split("\n") if line]
    os.makedirs(DATA_DIR, exist_ok=True)
    with open(cache, "w") as f:
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
