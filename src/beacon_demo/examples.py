"""Example requests, answered by the node itself. ``beacon-demo examples out.json`` writes them with their responses.

Each example asks one question at one or more access levels. Positions in requests are 0-based, as Beacon v2 uses
them: rs4244285 sits at VCF position 94781859 on chromosome 10, so the request says start=94781858.
"""

from __future__ import annotations

from urllib.parse import urlencode

TOKENS = {"PUBLIC": None, "REGISTERED": "registered-demo", "CONTROLLED": "controlled-demo"}
ALL_LEVELS = ("PUBLIC", "REGISTERED", "CONTROLLED")

EXAMPLES = [
    {"id": "info", "title": "What is this node?", "method": "GET", "path": "/info", "query": {},
     "levels": ("PUBLIC",),
     "why": "Every Beacon describes itself, so a network can list it."},
    {"id": "filtering-terms", "title": "What can I filter on?", "method": "GET", "path": "/filtering_terms",
     "query": {}, "levels": ("PUBLIC",),
     "why": "The phenotype and sex terms the node's people are described with (HPO and NCIT)."},
    {"id": "variant-yes-no", "title": "Is CYP2C19 rs4244285 (G>A) seen here?", "method": "GET", "path": "/g_variants",
     "query": {"referenceName": "10", "start": "94781858", "referenceBases": "G", "alternateBases": "A"},
     "levels": ("PUBLIC",),
     "why": "The basic Beacon question: a yes or no about one variant. No person-level data leaves the node."},
    {"id": "variant-frequency", "title": "How common is it in each dataset?", "method": "GET", "path": "/g_variants",
     "query": {"referenceName": "10", "start": "94781858", "referenceBases": "G", "alternateBases": "A",
               "requestedGranularity": "record"},
     "levels": ("PUBLIC",),
     "why": "An allele-frequency Beacon answers with counts and frequencies per dataset, not with people. "
            "These frequencies are made up."},
    {"id": "genotyped-not-seen", "title": "A variant that was looked for and not seen",
     "method": "GET", "path": "/g_variants",
     "query": {"referenceName": "13", "start": "48037781", "referenceBases": "AGGAGTC",
               "alternateBases": "AGGAGTCGGAGTC"},
     "levels": ("PUBLIC",),
     "why": "NUDT15 rs746071566, the duplication allele: the position was genotyped in both datasets and nobody "
            "carries it, so the answer is a real no."},
    {"id": "not-genotyped", "title": "An off-by-one position: nothing was genotyped there",
     "method": "GET", "path": "/g_variants",
     "query": {"referenceName": "10", "start": "94781859", "referenceBases": "G", "alternateBases": "A"},
     "levels": ("PUBLIC",),
     "why": "The same variant asked with its 1-based VCF position by mistake. A plain Beacon would say no; this "
            "node also says that nothing was genotyped there, so the no means nothing."},
    {"id": "deletion-vcf-style", "title": "One deletion, written as a VCF writes it", "method": "GET",
     "path": "/g_variants",
     "query": {"referenceName": "10", "start": "94949280", "referenceBases": "GA", "alternateBases": "G",
               "requestedGranularity": "count"},
     "levels": ("PUBLIC",),
     "why": "CYP2C9 rs9332131, a one-base deletion, with the base before it as an anchor."},
    {"id": "deletion-trimmed", "title": "The same deletion, written without the anchor base", "method": "GET",
     "path": "/g_variants",
     "query": {"referenceName": "10", "start": "94949281", "referenceBases": "A", "alternateBases": "-",
               "requestedGranularity": "count"},
     "levels": ("PUBLIC",),
     "why": "Both spellings find the same variant, because the node compares variants after trimming the bases "
            "both alleles share."},
    {"id": "gene", "title": "Which variants were seen in CYP2C9?", "method": "GET", "path": "/g_variants",
     "query": {"geneId": "CYP2C9", "requestedGranularity": "count"},
     "levels": ("PUBLIC",),
     "why": "A gene query. The count is of variants per dataset, not of people."},
    {"id": "people-with-phenotype", "title": "How many people have atrial fibrillation (HP:0005110)?",
     "method": "GET", "path": "/individuals",
     "query": {"filters": "HP:0005110", "requestedGranularity": "record", "limit": "2"},
     "levels": ALL_LEVELS,
     "why": "The same question at three access levels: the public gets yes or no, a registered user a count, "
            "an approved user the records."},
    {"id": "carriers-with-phenotype",
     "title": "How many people with hypertension (HP:0000822) carry CYP2C19 rs4244285?",
     "method": "POST", "path": "/g_variants/chr10-94781859-G-A/individuals",
     "body": {"meta": {"apiVersion": "v2.2.0"},
              "query": {"filters": [{"id": "HP:0000822"}], "requestedGranularity": "record",
                        "pagination": {"skip": 0, "limit": 2}}},
     "levels": ALL_LEVELS,
     "why": "A phenotype-linked question, sent as a POST. In this synthetic data, phenotypes were drawn "
            "independently of genotypes, so any number here is chance."},
    {"id": "small-count", "title": "Small counts are not given exactly",
     "method": "GET", "path": "/g_variants/chr16-31096368-C-T/individuals",
     "query": {"filters": "HP:0000822", "requestedGranularity": "count"},
     "levels": ("REGISTERED", "CONTROLLED"),
     "why": "VKORC1 rs9923231 carriers with hypertension. A count from 1 to 4 is reported as 5 below the "
            "controlled level, so a small group cannot be singled out."},
    {"id": "wrong-assembly", "title": "A question about another genome build is refused",
     "method": "GET", "path": "/g_variants",
     "query": {"referenceName": "10", "start": "94781858", "referenceBases": "G", "alternateBases": "A",
               "assemblyId": "GRCh37"},
     "levels": ("PUBLIC",),
     "why": "GRCh37 positions differ from GRCh38 ones. Answering no would be wrong, so the node returns an error."},
    {"id": "unknown-filter", "title": "An unknown filter is refused, not ignored",
     "method": "GET", "path": "/individuals", "query": {"filters": "HP:0001627"},
     "levels": ("PUBLIC",),
     "why": "Ignoring a filter would answer a broader question than the one asked."},
]


def run(beacon):
    out = []
    for example in EXAMPLES:
        query_string = urlencode(example.get("query", {}))
        for level in example["levels"]:
            token = TOKENS[level]
            headers = {"Authorization": f"Bearer {token}"} if token else {}
            status, response = beacon.handle(example["method"], example["path"], query_string,
                                             example.get("body"), headers)
            url = f"{beacon.base_url}{example['path']}" + (f"?{query_string}" if query_string else "")
            item = {"id": example["id"], "title": example["title"], "why": example["why"],
                    "accessLevel": level, "token": token, "method": example["method"], "url": url,
                    "status": status, "response": response}
            if "body" in example:
                item["body"] = example["body"]
            out.append(item)
    return out
