"""Everything that talks to NCBI: fetching the protein family the whole TP
works with, plus the local cache so the same records are not downloaded twice.

The family is beta-globin, one entry per species, and the ten entries are
named explicitly below rather than searched for. Two reasons for that choice
over "the top hit of a query per species", which is what this module used to
do and which quietly downloaded haptoglobin for the rabbit and a transcription
factor for the mouse:

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
(((((human, chimp), macaque), ((mouse, rat), rabbit)), cow), chicken, frog,
fish) -- primates and rodents as tight clades, fish at the root.
"""

import os
import json

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

# A globin is ~146 residues. The guard is not paranoia: the previous version
# of this module returned entries of 373, 347 and 969 aa and nothing noticed
# until the BLAST output looked strange. Anything outside this range means
# the accession is not what it is supposed to be, and that should stop the
# run rather than flow into a distance matrix.
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
