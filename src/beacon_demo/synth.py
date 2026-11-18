"""A synthetic cohort: made-up people with made-up genotypes at real GRCh38 positions.

Nothing here describes a real person or a real population. The positions and alleles are real (see positions.tsv);
the people, their sex, their phenotypes, their genotypes and the allele frequencies used to draw them are generated
from a fixed seed, so the same seed always gives the same files. Phenotypes are drawn independently of genotypes,
so the data holds no genotype-phenotype association of any kind.
"""

from __future__ import annotations

import math
import os
import random
from dataclasses import dataclass
from importlib import resources

from . import __version__

DEFAULT_SEED = 20261118
DATASETS = (
    ("synthetic-cohort-1", "Synthetic cohort 1", 200),
    ("synthetic-cohort-2", "Synthetic cohort 2", 100),
)
# (HPO term, label, share of made-up people given it). Labels checked against the HPO on 18 Nov 2026.
PHENOTYPES = (
    ("HP:0000822", "Hypertension", 0.25),
    ("HP:0005978", "Type II diabetes mellitus", 0.10),
    ("HP:0001513", "Obesity", 0.18),
    ("HP:0002099", "Asthma", 0.08),
    ("HP:0002076", "Migraine", 0.12),
    ("HP:0003124", "Hypercholesterolemia", 0.15),
    ("HP:0005110", "Atrial fibrillation", 0.04),
    ("HP:0001250", "Seizure", 0.02),
)
MISSING_RATE = 0.01  # share of genotypes written as ./.
CONTIG_LENGTHS = {  # GRCh38
    "chr1": 248956422, "chr4": 190214555, "chr6": 170805979, "chr7": 159345973, "chr10": 133797422,
    "chr12": 133275309, "chr13": 114364328, "chr16": 90338345, "chr19": 58617616,
}


@dataclass(frozen=True)
class Site:
    chrom: str
    pos: int
    id: str
    ref: str
    alts: tuple
    gene: str


@dataclass(frozen=True)
class Person:
    individual_id: str
    biosample_id: str
    dataset_id: str
    sex: str
    phenotypes: tuple


@dataclass
class Cohort:
    sites: list
    people: list
    genotypes: dict  # biosample id -> list of (allele, allele) tuples or None, one per site
    frequencies: list  # per site, the made-up frequency of each alternate allele used to draw genotypes


def load_sites():
    text = resources.files("beacon_demo").joinpath("positions.tsv").read_text(encoding="utf-8")
    sites = []
    for line in text.splitlines():
        if not line or line.startswith("#"):
            continue
        chrom, pos, vid, ref, alt, gene = line.split("\t")
        sites.append(Site(chrom, int(pos), vid, ref, tuple(alt.split(",")), gene))
    return sites


def _draw_allele(rng, freqs):
    r = rng.random()
    total = 0.0
    for i, f in enumerate(freqs, start=1):
        total += f
        if r < total:
            return i
    return 0


def make_cohort(seed=DEFAULT_SEED, datasets=DATASETS):
    rng = random.Random(seed)
    sites = load_sites()
    frequencies = []
    for site in sites:
        # One made-up frequency per alternate allele, spread on a log scale from 0.4% to 35%; extra alleles are rarer.
        first = math.exp(rng.uniform(math.log(0.004), math.log(0.35)))
        frequencies.append([round(first, 4)] + [round(first / 5, 4) for _ in site.alts[1:]])
    people = []
    n = 0
    for dataset_id, _name, size in datasets:
        for _ in range(size):
            n += 1
            sex = rng.choice(("female", "male"))
            phenotypes = tuple(term for term, _label, share in PHENOTYPES if rng.random() < share)
            people.append(Person(f"IND-{n:04d}", f"SMP-{n:04d}", dataset_id, sex, phenotypes))
    genotypes = {}
    for person in people:
        calls = []
        for freqs in frequencies:
            if rng.random() < MISSING_RATE:
                calls.append(None)
            else:
                calls.append(tuple(sorted((_draw_allele(rng, freqs), _draw_allele(rng, freqs)))))
        genotypes[person.biosample_id] = calls
    return Cohort(sites, people, genotypes, frequencies)


def _gt(call):
    return "./." if call is None else f"{call[0]}/{call[1]}"


def vcf_text(cohort, dataset_id):
    samples = [p.biosample_id for p in cohort.people if p.dataset_id == dataset_id]
    lines = [
        "##fileformat=VCFv4.2",
        f"##source=beacon-v2-synthetic-demo {__version__} (synthetic data: no real people)",
        "##reference=GRCh38",
    ]
    for chrom in sorted({s.chrom for s in cohort.sites}, key=lambda c: int(c[3:])):
        lines.append(f"##contig=<ID={chrom},length={CONTIG_LENGTHS[chrom]},assembly=GRCh38>")
    lines += [
        '##INFO=<ID=GENE,Number=1,Type=String,Description="Gene symbol">',
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
        "\t".join(["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"] + samples),
    ]
    for i, site in enumerate(cohort.sites):
        row = [site.chrom, str(site.pos), site.id, site.ref, ",".join(site.alts), ".", "PASS", f"GENE={site.gene}", "GT"]
        row += [_gt(cohort.genotypes[s][i]) for s in samples]
        lines.append("\t".join(row))
    return "\n".join(lines) + "\n"


def write(cohort, outdir, datasets=DATASETS):
    """Write the files the server reads: datasets.tsv, samples.tsv, terms.tsv and one VCF per dataset."""
    os.makedirs(outdir, exist_ok=True)
    written = []

    def put(name, text):
        path = os.path.join(outdir, name)
        with open(path, "w", encoding="utf-8", newline="\n") as fh:
            fh.write(text)
        written.append(path)

    rows = ["#id\tname\tvcf\tdescription"]
    for dataset_id, name, size in datasets:
        put(f"{dataset_id}.vcf", vcf_text(cohort, dataset_id))
        rows.append(f"{dataset_id}\t{name}\t{dataset_id}.vcf\t{size} made-up people with made-up genotypes at "
                    f"{len(cohort.sites)} real GRCh38 positions. Synthetic data: no real people.")
    put("datasets.tsv", "\n".join(rows) + "\n")
    rows = ["#individual_id\tbiosample_id\tdataset_id\tsex\tphenotypes"]
    for p in cohort.people:
        rows.append("\t".join([p.individual_id, p.biosample_id, p.dataset_id, p.sex, ";".join(p.phenotypes)]))
    put("samples.tsv", "\n".join(rows) + "\n")
    put("terms.tsv", "#id\tlabel\n" + "".join(f"{t}\t{label}\n" for t, label, _ in PHENOTYPES))
    return written
