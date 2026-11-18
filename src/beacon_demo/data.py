"""Load a data directory (or a synthetic cohort) into the in-memory index the Beacon answers from.

A data directory holds:

- ``datasets.tsv``: id, name, VCF file name, description; one row per dataset.
- ``samples.tsv``: individual id, biosample id (the VCF column name), dataset id, sex (female, male or empty),
  phenotypes (HPO terms separated by ``;``).
- ``terms.tsv``: id and label of every phenotype term used, for the filtering-terms endpoint.
- one VCF per dataset (plain or gzipped), GRCh38, with a GT field.

Everything stays in memory: this is a demonstration, not a database.
"""

from __future__ import annotations

import gzip
import os
import re
from dataclasses import dataclass, field

from . import synth
from .alleles import minimal, variant_type

_BASES = re.compile(r"^[ACGTN]+$")


@dataclass(frozen=True)
class Variant:
    vid: str
    chrom: str
    pos: int  # VCF position, 1-based
    ref: str
    alt: str
    rsid: str
    gene: str
    start: int  # 0-based interbase start of the VCF-style REF span
    end: int
    min_start: int  # minimal representation: shared bases trimmed, 0-based interbase
    min_ref: str
    min_alt: str
    vtype: str


@dataclass
class Stats:
    allele_number: int = 0
    allele_count: int = 0
    homozygous: int = 0
    heterozygous: int = 0
    hemizygous: int = 0
    carriers: dict = field(default_factory=dict)  # biosample id -> (copies of the allele, ploidy of the call)


@dataclass
class Dataset:
    id: str
    name: str
    description: str
    vcf: str
    biosamples: list = field(default_factory=list)
    stats: dict = field(default_factory=dict)  # variant id -> Stats
    sites: list = field(default_factory=list)  # (chrom, start, end, gene) of every VCF row: what was tested


@dataclass(frozen=True)
class Individual:
    id: str
    sex: str
    phenotypes: tuple
    biosample_id: str
    dataset_id: str


@dataclass
class Data:
    datasets: dict  # id -> Dataset, in file order
    variants: dict  # variant id -> Variant
    individuals: dict  # id -> Individual
    biosample_owner: dict  # biosample id -> individual id
    terms: dict  # phenotype term -> label
    skipped: list = field(default_factory=list)  # notes about input rows that were not loaded

    def variants_of(self, biosample_id):
        dataset = self.datasets[self.individuals[self.biosample_owner[biosample_id]].dataset_id]
        return [vid for vid, st in dataset.stats.items() if biosample_id in st.carriers]


def variant_id(chrom, pos, ref, alt):
    return f"{chrom}-{pos}-{ref}-{alt}"


def _parse_gt(text):
    gt = text.split(":", 1)[0].strip()
    parts = re.split(r"[/|]", gt)
    if any(p in ("", ".") for p in parts):
        return None
    try:
        return tuple(int(p) for p in parts)
    except ValueError:
        return None


def _rows(path):
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt", encoding="utf-8") as fh:
        for line in fh:
            yield line.rstrip("\n")


def _tsv(path):
    out = []
    for line in _rows(path):
        if line and not line.startswith("#"):
            out.append(line.split("\t"))
    return out


def _add_vcf(data, dataset, lines):
    samples = None
    for line in lines:
        if line.startswith("##") or not line:
            continue
        cols = line.split("\t")
        if line.startswith("#CHROM"):
            samples = cols[9:]
            missing = [s for s in samples if s not in data.biosample_owner]
            if missing:
                raise ValueError(f"{dataset.vcf}: VCF columns not in samples.tsv: {', '.join(missing[:5])}")
            dataset.biosamples = samples
            continue
        if samples is None:
            raise ValueError(f"{dataset.vcf}: data row before the #CHROM header")
        chrom, pos, rsid, ref, alts = cols[0], int(cols[1]), cols[2], cols[3].upper(), cols[4].upper().split(",")
        chrom = chrom if chrom.startswith("chr") else "chr" + chrom
        gene = ""
        for item in cols[7].split(";"):
            if item.startswith("GENE="):
                gene = item[5:]
        if not _BASES.match(ref):
            data.skipped.append(f"{dataset.id}: {chrom}:{pos} REF {ref!r} is not plain bases")
            continue
        fmt = cols[8].split(":")
        if fmt[0] != "GT":
            data.skipped.append(f"{dataset.id}: {chrom}:{pos} has no GT as its first FORMAT field")
            continue
        calls = [_parse_gt(c) for c in cols[9:]]
        dataset.sites.append((chrom, pos - 1, pos - 1 + len(ref), gene))
        for index, alt in enumerate(alts, start=1):
            if not _BASES.match(alt):
                data.skipped.append(f"{dataset.id}: {chrom}:{pos} ALT {alt!r} is not plain bases")
                continue
            vid = variant_id(chrom, pos, ref, alt)
            if vid not in data.variants:
                m_start, m_ref, m_alt = minimal(pos - 1, ref, alt)
                data.variants[vid] = Variant(vid, chrom, pos, ref, alt, "" if rsid == "." else rsid, gene,
                                             pos - 1, pos - 1 + len(ref), m_start, m_ref, m_alt,
                                             variant_type(m_ref, m_alt))
            st = dataset.stats.setdefault(vid, Stats())
            for biosample, call in zip(samples, calls):
                if call is None:
                    continue
                st.allele_number += len(call)
                copies = sum(1 for a in call if a == index)
                if copies:
                    st.allele_count += copies
                    st.carriers[biosample] = (copies, len(call))
                    if len(call) == 1:
                        st.hemizygous += 1
                    elif len(call) == 2:
                        if copies == 2:
                            st.homozygous += 1
                        else:
                            st.heterozygous += 1


def _finish(datasets, individuals, terms, vcf_lines):
    data = Data({d.id: d for d in datasets}, {}, {i.id: i for i in individuals},
                {i.biosample_id: i.id for i in individuals}, terms)
    for d in datasets:
        _add_vcf(data, d, vcf_lines(d))
    return data


def from_cohort(cohort):
    """Index a synthetic cohort without writing files."""
    datasets = [Dataset(i, n, f"{size} made-up people with made-up genotypes at {len(cohort.sites)} real GRCh38 "
                        f"positions. Synthetic data: no real people.", f"{i}.vcf")
                for i, n, size in synth.DATASETS]
    individuals = [Individual(p.individual_id, p.sex, p.phenotypes, p.biosample_id, p.dataset_id) for p in cohort.people]
    terms = {t: label for t, label, _ in synth.PHENOTYPES}
    return _finish(datasets, individuals, terms, lambda d: synth.vcf_text(cohort, d.id).splitlines())


def load(directory):
    """Index a data directory written by ``beacon-demo synth`` or laid out the same way."""
    datasets = [Dataset(r[0], r[1], r[3] if len(r) > 3 else "", r[2]) for r in _tsv(os.path.join(directory, "datasets.tsv"))]
    known = {d.id for d in datasets}
    individuals = []
    for r in _tsv(os.path.join(directory, "samples.tsv")):
        r += [""] * (5 - len(r))
        if r[2] not in known:
            raise ValueError(f"samples.tsv: {r[0]} names an unknown dataset {r[2]!r}")
        sex = r[3].strip().lower()
        if sex not in ("female", "male", ""):
            raise ValueError(f"samples.tsv: {r[0]} has sex {r[3]!r}; use female, male or leave it empty")
        individuals.append(Individual(r[0], sex, tuple(t for t in r[4].split(";") if t), r[1], r[2]))
    terms = {r[0]: r[1] for r in _tsv(os.path.join(directory, "terms.tsv"))}
    unknown = sorted({t for i in individuals for t in i.phenotypes} - set(terms))
    if unknown:
        raise ValueError(f"samples.tsv uses terms missing from terms.tsv: {', '.join(unknown)}")
    return _finish(datasets, individuals, terms, lambda d: _rows(os.path.join(directory, d.vcf)))
