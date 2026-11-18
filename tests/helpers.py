import functools

from beacon_demo import data, synth
from beacon_demo.beacon import Beacon


@functools.lru_cache(maxsize=None)
def cohort():
    return synth.make_cohort()


def beacon():
    return Beacon(data.from_cohort(cohort()))


def ask(b, path, query="", level="PUBLIC", method="GET", body=None):
    token = {"PUBLIC": None, "REGISTERED": "registered-demo", "CONTROLLED": "controlled-demo"}[level]
    headers = {"Authorization": f"Bearer {token}"} if token else {}
    return b.handle(method, path, query, body, headers)


def vcf_counts(dataset_id, chrom, pos, alt):
    """Allele count, allele number and carriers, read straight from the synthetic VCF text."""
    lines = synth.vcf_text(cohort(), dataset_id).splitlines()
    samples = next(line for line in lines if line.startswith("#CHROM")).split("\t")[9:]
    for line in lines:
        cols = line.split("\t")
        if line.startswith("#") or cols[0] != chrom or int(cols[1]) != pos:
            continue
        index = cols[4].split(",").index(alt) + 1
        ac = an = 0
        carriers = {}
        for sample, gt in zip(samples, cols[9:]):
            if "." in gt:
                continue
            alleles = [int(a) for a in gt.split("/")]
            an += len(alleles)
            copies = alleles.count(index)
            ac += copies
            if copies:
                carriers[sample] = copies
        return ac, an, carriers
    raise KeyError((chrom, pos, alt))
