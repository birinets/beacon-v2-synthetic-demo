"""Answer Beacon v2 requests from the in-memory index. The HTTP layer is in server.py.

Access levels (a demonstration of the idea, not security):

- PUBLIC (no token): yes/no answers about people and samples; variant records with aggregate frequencies only, as an
  allele-frequency Beacon gives them. A variant query combined with filters on people answers yes/no.
- REGISTERED (``Authorization: Bearer registered-demo``): counts. Counts from 1 to 4 are reported as 5.
- CONTROLLED (``Authorization: Bearer controlled-demo``): records, including which sample carries which variant.

The Beacon never answers a question it cannot read: unknown parameters, filters, assemblies and query types are
refused with an error instead of being ignored, and every variant query says whether the queried region was
genotyped in each dataset, because "exists": false at a position nobody was genotyped at says nothing about anyone.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from urllib.parse import parse_qs

from . import __version__, spec
from .alleles import clean, length, minimal

SMALL_COUNT = 5
DEMO_TOKENS = {"registered-demo": "REGISTERED", "controlled-demo": "CONTROLLED"}
NOTICE = "Synthetic data: made-up people and genotypes at real GRCh38 positions. No real people."
UNKNOWN_SEX = {"id": "NCIT:C17998", "label": "unknown"}

_VARIANT_PARAMS = ("referenceName", "start", "end", "referenceBases", "alternateBases", "variantType",
                   "variantMinLength", "variantMaxLength", "geneId", "assemblyId", "aminoacidChange",
                   "genomicAlleleShortForm", "mateName")
_GET_ALIASES = {"geneid": "geneId", "allele": "genomicAlleleShortForm"}
_GET_CONTROL = ("requestedGranularity", "skip", "limit", "includeResultsetResponses", "testMode", "filters",
                "requestedSchema")
_BASES = re.compile(r"^[ACGTN]*$")
_GRCH38 = re.compile(r"^(grch38|hg38)(\.p\d+)?$", re.IGNORECASE)


class BeaconError(Exception):
    def __init__(self, status, message):
        super().__init__(message)
        self.status = status
        self.message = message


@dataclass
class Request:
    api_version: str = spec.API_VERSION
    requested_schemas: list = field(default_factory=list)
    granularity: str = "boolean"
    include: str = "HIT"
    skip: int = 0
    limit: int = spec.DEFAULT_LIMIT
    test_mode: bool = False
    params: dict = field(default_factory=dict)
    filters: list = field(default_factory=list)  # dicts with at least "id"


# ---------------------------------------------------------------- request parsing

def _int(name, value, minimum=0):
    try:
        number = int(value)
    except (TypeError, ValueError):
        raise BeaconError(400, f"{name} must be an integer, not {value!r}")
    if number < minimum:
        raise BeaconError(400, f"{name} must be at least {minimum}")
    return number


def _int_list(name, value):
    items = value if isinstance(value, list) else [value] if isinstance(value, int) else str(value).split(",")
    numbers = [_int(name, v) for v in items]
    if not 1 <= len(numbers) <= 2:
        raise BeaconError(400, f"{name} takes one or two values")
    return numbers


def _bool(name, value):
    if isinstance(value, bool):
        return value
    if str(value).lower() in ("true", "1"):
        return True
    if str(value).lower() in ("false", "0"):
        return False
    raise BeaconError(400, f"{name} must be true or false")


def _check_control(req):
    if req.granularity not in spec.GRANULARITIES:
        raise BeaconError(400, f"requestedGranularity must be one of {', '.join(spec.GRANULARITIES)}")
    if req.include not in spec.INCLUDE_RESULTSET_RESPONSES:
        raise BeaconError(400, f"includeResultsetResponses must be one of {', '.join(spec.INCLUDE_RESULTSET_RESPONSES)}")


def _variant_param(req, name, value):
    if name in ("start", "end"):
        req.params[name] = _int_list(name, value)
    elif name in ("variantMinLength", "variantMaxLength"):
        req.params[name] = _int(name, value)
    else:
        req.params[name] = value if isinstance(value, str) else str(value)


def parse_get(query_string):
    req = Request()
    for key, values in parse_qs(query_string, keep_blank_values=True).items():
        value = values[-1]
        name = _GET_ALIASES.get(key, key)
        if name in _VARIANT_PARAMS:
            _variant_param(req, name, value)
        elif name == "requestedGranularity":
            req.granularity = value
        elif name == "skip":
            req.skip = _int("skip", value)
        elif name == "limit":
            req.limit = _int("limit", value)
        elif name == "includeResultsetResponses":
            req.include = value
        elif name == "testMode":
            req.test_mode = _bool("testMode", value)
        elif name == "filters":
            req.filters = [{"id": f.strip()} for v in values for f in v.split(",") if f.strip()]
        elif name == "requestedSchema":
            req.requested_schemas = [{"schema": value}]
        else:
            raise BeaconError(400, f"unknown query parameter {key!r}; supported: "
                                   f"{', '.join(_VARIANT_PARAMS[:-3] + _GET_CONTROL)}")
    _check_control(req)
    return req


def parse_post(body):
    if body is None:
        return Request()
    if not isinstance(body, dict):
        raise BeaconError(400, "the request body must be a JSON object")
    meta = body.get("meta") or {}
    query = body.get("query") or {}
    if not isinstance(meta, dict) or not isinstance(query, dict):
        raise BeaconError(400, "meta and query must be JSON objects")
    req = Request()
    req.api_version = str(meta.get("apiVersion", spec.API_VERSION))
    req.requested_schemas = meta.get("requestedSchemas") or []
    req.granularity = query.get("requestedGranularity", meta.get("requestedGranularity", "boolean"))
    pagination = query.get("pagination") or {}
    req.skip = _int("pagination.skip", pagination.get("skip", 0))
    req.limit = _int("pagination.limit", pagination.get("limit", spec.DEFAULT_LIMIT))
    req.include = query.get("includeResultsetResponses", "HIT")
    req.test_mode = _bool("testMode", query.get("testMode", False))
    params = query.get("requestParameters") or {}
    if not isinstance(params, dict):
        raise BeaconError(400, "requestParameters must be one JSON object; lists of variants are not supported here")
    params = {k: v for k, v in params.items() if k != "$schema"}
    sets = [k for k, v in params.items() if isinstance(v, dict)]
    if sets:  # the framework's form, {"g_variant": {...}}; the flat form is accepted too
        if sets != ["g_variant"] or len(params) > 1:
            raise BeaconError(400, "requestParameters: give the variant parameters either under g_variant or flat; "
                                   f"this node reads no other parameter set ({', '.join(sets)} given)")
        params = params["g_variant"]
    for name, value in params.items():
        if name not in _VARIANT_PARAMS:
            raise BeaconError(400, f"unknown request parameter {name!r}")
        _variant_param(req, name, value)
    filters = query.get("filters") or []
    if not isinstance(filters, list):
        raise BeaconError(400, "filters must be a list")
    for f in filters:
        if isinstance(f, str):
            req.filters.append({"id": f})
        elif isinstance(f, dict) and isinstance(f.get("id"), str):
            req.filters.append(dict(f))
        else:
            raise BeaconError(400, "each filter must be a term id or an object with an id")
    _check_control(req)
    return req


def access_level(headers):
    auth = headers.get("authorization", "").strip()
    if not auth:
        return "PUBLIC"
    scheme, _, token = auth.partition(" ")
    level = DEMO_TOKENS.get(token.strip()) if scheme.lower() == "bearer" else None
    if level is None:
        raise BeaconError(401, "unknown access token; this demo knows only registered-demo and controlled-demo")
    return level


# ---------------------------------------------------------------- variant queries

@dataclass
class Query:
    kind: str  # all, gene, sequence, range, bracket, single
    match: object  # Variant -> bool
    tested: object  # Dataset -> bool: does the dataset hold genotypes in the queried region?


def _chrom(name):
    text = str(name).strip()
    for chrom, accession in spec.REFSEQ_GRCH38.items():
        if text == accession:
            return chrom
    bare = text[3:] if text.lower().startswith("chr") else text
    if bare.isdigit() and 1 <= int(bare) <= 22:
        return f"chr{int(bare)}"
    if bare.upper() in ("X", "Y"):
        return f"chr{bare.upper()}"
    raise BeaconError(400, f"referenceName {name!r} is not a GRCh38 chromosome; use 1-22, X, Y, chr10 or NC_000010.11")


def _overlaps(dataset, chrom, start, end):
    return any(c == chrom and s < end and e > start for c, s, e, _ in dataset.sites)


def compile_query(params):
    p = dict(params)
    for name in ("aminoacidChange", "genomicAlleleShortForm", "mateName"):
        if name in p:
            raise BeaconError(400, f"{name} queries are not supported by this demo")
    assembly = p.get("assemblyId")
    if assembly is not None and not _GRCH38.match(assembly):
        raise BeaconError(400, f"assemblyId {assembly!r}: this Beacon holds GRCh38 data only")
    vtype = p.get("variantType")
    min_len, max_len = p.get("variantMinLength"), p.get("variantMaxLength")
    symbol = p["geneId"].strip().upper() if p.get("geneId") is not None else None

    def extra(v):
        if symbol is not None and v.gene.upper() != symbol:
            return False
        if vtype is not None:
            wanted = vtype.upper()
            if not (v.vtype == wanted or (wanted == "INDEL" and v.vtype in ("DEL", "INS", "INDEL"))):
                return False
        size = length(v.min_ref, v.min_alt)
        return (min_len is None or size >= min_len) and (max_len is None or size <= max_len)

    gene, start, end, name = p.get("geneId"), p.get("start"), p.get("end"), p.get("referenceName")
    has_alleles = "referenceBases" in p or "alternateBases" in p
    if gene is None and start is None and end is None and name is None:
        if has_alleles:
            raise BeaconError(400, "referenceBases and alternateBases need referenceName and start")
        return Query("all", extra, None)
    if gene is not None and start is None:  # a gene query; with start (and end) it narrows a position query instead
        chrom = _chrom(name) if name is not None else None
        return Query("gene", lambda v: (chrom is None or v.chrom == chrom) and extra(v),
                     lambda d: any(site[3].upper() == symbol for site in d.sites))
    if name is None or start is None:
        raise BeaconError(400, "a sequence, range or bracket query needs referenceName and start (0-based)")
    chrom = _chrom(name)
    ends = end or []
    if len(start) == 2 or len(ends) == 2:
        if len(start) != 2 or len(ends) != 2:
            raise BeaconError(400, "a bracket query needs two start values and two end values")
        s1, s2 = start
        e1, e2 = ends
        if s1 > s2 or e1 > e2:
            raise BeaconError(400, "bracket values must be given as minimum,maximum")
        return Query("bracket",
                     lambda v: v.chrom == chrom and s1 <= v.start <= s2 and e1 <= v.end <= e2 and extra(v),
                     lambda d: _overlaps(d, chrom, s1, e2))
    s = start[0]
    alt = clean(p.get("alternateBases"))
    ref = clean(p.get("referenceBases"))
    for label, bases in (("referenceBases", ref), ("alternateBases", alt)):
        if bases is not None and not _BASES.match(bases):
            raise BeaconError(400, f"{label}: only A, C, G, T, N, '-' and '.' are supported by this demo")
    if ends:
        e = ends[0]
        if e <= s:
            raise BeaconError(400, "end must be greater than start")

        def in_range(v):
            if not (v.chrom == chrom and v.start < e and v.end > s and extra(v)):
                return False
            if alt is None:
                return True
            if alt == "N":
                return v.vtype == "SNP"
            return alt in (v.alt, v.min_alt)
        return Query("range", in_range, lambda d: _overlaps(d, chrom, s, e))
    if alt is None:
        raise BeaconError(400, "a sequence query needs alternateBases and referenceBases; give end for a range query")
    if alt == "N":
        if ref is not None and len(ref) > 1:
            raise BeaconError(400, "alternateBases N stands for any single base, so referenceBases must be one base")
        match = (lambda v: v.chrom == chrom and v.vtype == "SNP" and v.min_start == s
                 and ref in (None, "N", v.min_ref) and extra(v))
        span = (s, s + 1)
    elif ref is None:
        match = lambda v: v.chrom == chrom and v.start == s and v.alt == alt and extra(v)
        span = (s, s + 1)
    else:
        if "N" in ref or "N" in alt:
            raise BeaconError(400, "N is accepted only as the whole alternateBases of a single-base query")
        key = minimal(s, ref, alt)
        if not key[1] and not key[2]:
            raise BeaconError(400, "referenceBases and alternateBases are the same: that is not a variant")
        match = lambda v: v.chrom == chrom and (v.min_start, v.min_ref, v.min_alt) == key and extra(v)
        span = (s, s + max(len(ref), 1))
    return Query("sequence", match, lambda d: _overlaps(d, chrom, *span))


# ---------------------------------------------------------------- the Beacon

def _lower(requested, ceiling):
    return spec.GRANULARITIES[min(spec.GRANULARITIES.index(requested), spec.GRANULARITIES.index(ceiling))]


def max_granularity(level, entry_type, person_level):
    """The most detail an access level may see. person_level: the answer depends on who carries what."""
    if entry_type == "dataset":
        return "record"
    if entry_type == "genomicVariant" and not person_level:
        return "record"  # aggregate frequencies only below CONTROLLED, as in an allele-frequency Beacon
    return {"PUBLIC": "boolean", "REGISTERED": "count", "CONTROLLED": "record"}[level]


class Beacon:
    def __init__(self, data, base_url="http://127.0.0.1:5050/api", beacon_id="org.example.beacon.synthetic-demo",
                 name="Beacon v2 demo node (synthetic data)", organization=None, data_version="2026-11-18",
                 created="2026-11-18T00:00:00+00:00"):
        self.data = data
        self.base_url = base_url.rstrip("/")
        self.beacon_id = beacon_id
        self.name = name
        self.organization = organization or {
            "id": "org.example", "name": "Example organisation (replace with the node's operator)",
            "welcomeUrl": "https://example.org"}
        self.data_version = data_version
        self.created = created

    # ---- entry point
    def error(self, status, message):
        """An error response for problems found before a request reaches handle(), such as unreadable JSON."""
        return status, {"meta": self._meta(Request(), "boolean", []),
                        "error": {"errorCode": status, "errorMessage": message}}

    def handle(self, method, path, query_string="", body=None, headers=None):
        req = Request()
        try:
            headers = {k.lower(): v for k, v in (headers or {}).items()}
            if method == "GET":
                req = parse_get(query_string or "")
            elif method == "POST":
                req = parse_post(body)
            else:
                raise BeaconError(405, f"method {method} is not supported; use GET or POST")
            level = access_level(headers)
            return 200, self._route([p for p in path.split("/") if p], req, level)
        except BeaconError as err:
            return err.status, {
                "meta": self._meta(req, "boolean", []),
                "error": {"errorCode": err.status, "errorMessage": err.message},
            }

    def _route(self, parts, req, level):
        if not parts or parts == ["info"]:
            return self.info()
        head, rest = parts[0], parts[1:]
        informational = {"service-info": self.service_info, "configuration": self.configuration, "map": self.beacon_map,
                         "entry_types": self.entry_types, "filtering_terms": self.filtering_terms}
        if head in informational and not rest:
            return informational[head]()
        datasets = list(self.data.datasets.values())
        if head == "datasets":
            if not rest:
                return self._datasets(req, datasets)
            dataset = self.data.datasets.get(rest[0])
            if dataset is None:
                raise BeaconError(404, f"no dataset {rest[0]!r}")
            if len(rest) == 1:
                return self._datasets(req, [dataset])
            if len(rest) == 2 and rest[1] in ("g_variants", "individuals", "biosamples"):
                return self._collection(rest[1], req, level, [dataset])
        elif head in ("g_variants", "individuals", "biosamples"):
            if not rest:
                return self._collection(head, req, level, datasets)
            if len(rest) == 1:
                return self._single(head, rest[0], req, level)
            if len(rest) == 2 and rest[1] in ("g_variants", "individuals", "biosamples") and rest[1] != head:
                return self._lookup(head, rest[0], rest[1], req, level)
        raise BeaconError(404, f"no endpoint /{'/'.join(parts)}; see {self.base_url}/map")

    # ---- filters on people
    def _people(self, filters):
        """Individual ids matching every filter, or None when there are no filters."""
        if not filters:
            return None
        chosen = set(self.data.individuals)
        sex_by_term = {t["id"]: sex for sex, t in spec.SEX_TERMS.items()}
        for f in filters:
            if "operator" in f or "value" in f:
                raise BeaconError(400, f"filter {f['id']!r}: alphanumeric filters are not supported by this demo")
            scope = f.get("scope")
            if scope is not None and scope not in ("individual", "individuals"):
                raise BeaconError(400, f"filter {f['id']!r}: this demo applies filters to individuals only")
            if f.get("similarity", "exact") != "exact":
                raise BeaconError(400, f"filter {f['id']!r}: only exact matching is supported (no similarity)")
            term = f["id"]
            if term in sex_by_term:
                chosen = {i for i in chosen if self.data.individuals[i].sex == sex_by_term[term]}
            elif term in self.data.terms:
                chosen = {i for i in chosen if term in self.data.individuals[i].phenotypes}
            else:
                raise BeaconError(400, f"unknown filter {term!r}; see {self.base_url}/filtering_terms")
        return chosen

    def _biosamples_of(self, individual_ids):
        return None if individual_ids is None else {self.data.individuals[i].biosample_id for i in individual_ids}

    # ---- collections and lookups
    def _variant_hits(self, query, datasets, samples):
        hits = {}
        for d in datasets:
            ids = []
            for vid, st in d.stats.items():
                if not st.carriers or not query.match(self.data.variants[vid]):
                    continue
                if samples is not None and not any(b in samples for b in st.carriers):
                    continue
                ids.append(vid)
            hits[d.id] = ids
        return hits

    def _coverage(self, query, datasets):
        if query.tested is None:
            return None
        return {d.id: bool(query.tested(d)) for d in datasets}

    def _collection(self, endpoint, req, level, datasets):
        people = self._people(req.filters)
        samples = self._biosamples_of(people)
        if endpoint == "g_variants":
            query = compile_query(req.params)
            hits = self._variant_hits(query, datasets, samples)
            return self._answer(req, level, "genomicVariant", hits, samples, people is not None,
                                self._coverage(query, datasets))
        coverage = None
        if req.params:  # people or samples who carry a variant matching the request parameters
            query = compile_query(req.params)
            variant_hits = self._variant_hits(query, datasets, samples)
            carriers = set()
            for d in datasets:
                for vid in variant_hits[d.id]:
                    carriers.update(d.stats[vid].carriers)
            samples = carriers if samples is None else samples & carriers
            coverage = self._coverage(query, datasets)
        hits = {}
        for d in datasets:
            chosen = [b for b in d.biosamples if samples is None or b in samples]
            hits[d.id] = chosen if endpoint == "biosamples" else [self.data.biosample_owner[b] for b in chosen]
        entry_type = "biosample" if endpoint == "biosamples" else "individual"
        return self._answer(req, level, entry_type, hits, samples, True, coverage)

    def _single(self, endpoint, entry_id, req, level):
        people = self._people(req.filters)
        samples = self._biosamples_of(people)
        datasets = list(self.data.datasets.values())
        if endpoint == "g_variants":
            if req.params:
                raise BeaconError(400, "request parameters do not apply to a single variant")
            hits = {d.id: [entry_id] if entry_id in d.stats and d.stats[entry_id].carriers and (
                samples is None or any(b in samples for b in d.stats[entry_id].carriers)) else [] for d in datasets}
            return self._answer(req, level, "genomicVariant", hits, samples, people is not None, None)
        if endpoint == "individuals":
            person = self.data.individuals.get(entry_id)
            ok = person is not None and (people is None or entry_id in people)
            hits = {d.id: [entry_id] if ok and person.dataset_id == d.id else [] for d in datasets}
            return self._answer(req, level, "individual", hits, samples, True, None)
        owner = self.data.biosample_owner.get(entry_id)
        ok = owner is not None and (people is None or owner in people)
        hits = {d.id: [entry_id] if ok and self.data.individuals[owner].dataset_id == d.id else [] for d in datasets}
        return self._answer(req, level, "biosample", hits, samples, True, None)

    def _lookup(self, endpoint, entry_id, target, req, level):
        datasets = list(self.data.datasets.values())
        people = self._people(req.filters)
        if endpoint == "g_variants":  # people or samples carrying one variant
            if req.params:
                raise BeaconError(400, "request parameters do not apply to a single variant")
            samples = self._biosamples_of(people)
            hits = {}
            for d in datasets:
                st = d.stats.get(entry_id)
                carriers = [b for b in d.biosamples if st and b in st.carriers and (samples is None or b in samples)]
                hits[d.id] = carriers if target == "biosamples" else [self.data.biosample_owner[b] for b in carriers]
            entry_type = "biosample" if target == "biosamples" else "individual"
            return self._answer(req, level, entry_type, hits, samples, True, None)
        if endpoint == "individuals":
            person = self.data.individuals.get(entry_id)
            biosample = person.biosample_id if person and (people is None or entry_id in people) else None
        else:
            owner = self.data.biosample_owner.get(entry_id)
            biosample = entry_id if owner and (people is None or owner in people) else None
        if target == "biosamples":  # the sample of one individual
            if req.params:
                raise BeaconError(400, "request parameters do not apply to this lookup")
            hits = {d.id: [biosample] if biosample in d.biosamples else [] for d in datasets}
            return self._answer(req, level, "biosample", hits, {biosample} if biosample else set(), True, None)
        if target == "individuals":  # the individual a sample came from
            if req.params:
                raise BeaconError(400, "request parameters do not apply to this lookup")
            hits = {d.id: [self.data.biosample_owner[biosample]] if biosample in d.biosamples else [] for d in datasets}
            return self._answer(req, level, "individual", hits, {biosample} if biosample else set(), True, None)
        query = compile_query(req.params)  # the variants one person carries, optionally narrowed by a query
        samples = {biosample} if biosample else set()
        hits = self._variant_hits(query, datasets, samples)
        return self._answer(req, level, "genomicVariant", hits, samples, True, self._coverage(query, datasets))

    # ---- responses
    def _summary(self, req):
        summary = {
            "apiVersion": req.api_version,
            "requestedSchemas": req.requested_schemas,
            "pagination": {"skip": req.skip, "limit": req.limit},
            "requestedGranularity": req.granularity,
            "includeResultsetResponses": req.include,
            "testMode": req.test_mode,
        }
        if req.filters:
            summary["filters"] = [f["id"] for f in req.filters]
        if req.params:  # keyed by parameter set, as the framework's requestParameters schema requires
            summary["requestParameters"] = {"g_variant": req.params}
        return summary

    def _meta(self, req, granularity, schemas):
        return {
            "beaconId": self.beacon_id,
            "apiVersion": spec.API_VERSION,
            "returnedGranularity": granularity,
            "returnedSchemas": schemas,
            "receivedRequestSummary": self._summary(req),
            "testMode": req.test_mode,
        }

    def _info_meta(self, entity, schema):
        return {"beaconId": self.beacon_id, "apiVersion": spec.API_VERSION,
                "returnedSchemas": [{"entityType": entity, "schema": schema}]}

    def _answer(self, req, level, entry_type, hits, samples, person_level, coverage):
        ceiling = max_granularity(level, entry_type, person_level)
        granularity = _lower(req.granularity, ceiling)
        total = sum(len(ids) for ids in hits.values())
        out = {"meta": self._meta(req, granularity, [spec.returned_schema(entry_type)])}
        summary = {"exists": total > 0}
        if granularity == "count":
            sensitive = person_level and level != "CONTROLLED"
            if sensitive and 0 < total < SMALL_COUNT:
                summary.update(numTotalResults=SMALL_COUNT, countPrecision="rounded", countAdjustedTo=[SMALL_COUNT])
            else:
                summary.update(numTotalResults=total, countPrecision="exact")
        elif granularity == "record":
            summary["numTotalResults"] = total
        out["responseSummary"] = summary
        if granularity == "record":
            render = {"genomicVariant": lambda d, x: self._variant_record(d, x, level, samples),
                      "individual": lambda d, x: self._individual_record(x),
                      "biosample": lambda d, x: self._biosample_record(x)}[entry_type]
            sets = []
            for dataset_id, ids in hits.items():
                hit = bool(ids)
                if req.include == "NONE" or (req.include == "HIT" and not hit) or (req.include == "MISS" and hit):
                    continue
                page = ids[req.skip:req.skip + req.limit] if req.limit else ids[req.skip:]
                result_set = {"id": dataset_id, "setType": "dataset", "exists": hit, "resultsCount": len(ids),
                              "results": [render(dataset_id, x) for x in page]}
                if coverage is not None:
                    result_set["info"] = {"queriedRegionGenotyped": coverage[dataset_id]}
                sets.append(result_set)
            out["response"] = {"resultSets": sets}
        info = {"notice": NOTICE, "accessLevel": level}
        if granularity != req.granularity:
            info["granularityNote"] = (f"{req.granularity} was requested; {granularity} is the most this access "
                                       f"level gets for this question")
        if coverage is not None:
            info["queriedRegionGenotyped"] = coverage
            if not all(coverage.values()):
                info["coverageNote"] = ("A dataset marked false holds no genotypes in the queried region, so "
                                        "\"exists\": false says nothing about its people there.")
        out["info"] = info
        return out

    def _variant_record(self, dataset_id, vid, level, samples):
        v = self.data.variants[vid]
        st = self.data.datasets[dataset_id].stats[vid]
        accession = spec.REFSEQ_GRCH38[v.chrom]
        record = {
            "variantInternalId": v.vid,
            "variation": {
                "variantType": v.vtype,
                "referenceBases": v.ref,
                "alternateBases": v.alt,
                "location": {
                    "type": "SequenceLocation",
                    "sequence_id": f"refseq:{accession}",
                    "interval": {"type": "SequenceInterval",
                                 "start": {"type": "Number", "value": v.start},
                                 "end": {"type": "Number", "value": v.end}},
                },
            },
        }
        identifiers = {}
        if v.vtype == "SNP" and len(v.ref) == 1:
            identifiers["genomicHGVSId"] = f"{accession}:g.{v.pos}{v.ref}>{v.alt}"
        if v.rsid:
            identifiers["variantAlternativeIds"] = [
                {"id": f"dbSNP:{v.rsid}", "reference": f"https://www.ncbi.nlm.nih.gov/snp/{v.rsid}"}]
        if identifiers:
            record["identifiers"] = identifiers
        if v.gene:
            record["molecularAttributes"] = {"geneIds": [v.gene]}
        frequency = {
            "population": dataset_id,
            "alleleFrequency": round(st.allele_count / st.allele_number, 6) if st.allele_number else 0.0,
            "alleleCount": st.allele_count,
            "alleleNumber": st.allele_number,
            "alleleCountHomozygous": st.homozygous,
            "alleleCountHeterozygous": st.heterozygous,
        }
        if st.hemizygous:
            frequency["alleleCountHemizygous"] = st.hemizygous
        record["frequencyInPopulations"] = [{
            "source": "Synthetic genotypes from beacon-v2-synthetic-demo (not a real population)",
            "sourceReference": f"{self.base_url}/datasets/{dataset_id}",
            "version": self.data_version,
            "frequencies": [frequency],
        }]
        if level == "CONTROLLED":
            cases = []
            for biosample, (copies, ploidy) in st.carriers.items():
                if samples is not None and biosample not in samples:
                    continue
                case = {"biosampleId": biosample, "individualId": self.data.biosample_owner[biosample]}
                if (copies, ploidy) in spec.ZYGOSITY_TERMS:
                    case["zygosity"] = spec.ZYGOSITY_TERMS[(copies, ploidy)]
                cases.append(case)
            record["caseLevelData"] = cases
        return record

    def _individual_record(self, individual_id):
        person = self.data.individuals[individual_id]
        record = {"id": person.id, "sex": spec.SEX_TERMS.get(person.sex, UNKNOWN_SEX)}
        if person.phenotypes:
            record["phenotypicFeatures"] = [{"featureType": {"id": t, "label": self.data.terms[t]}}
                                            for t in person.phenotypes]
        record["info"] = {"datasetId": person.dataset_id, "biosampleId": person.biosample_id}
        return record

    def _biosample_record(self, biosample_id):
        owner = self.data.biosample_owner[biosample_id]
        return {"id": biosample_id, "individualId": owner, "biosampleStatus": spec.BIOSAMPLE_STATUS,
                "sampleOriginType": spec.SAMPLE_ORIGIN_TYPE,
                "info": {"datasetId": self.data.individuals[owner].dataset_id}}

    def _dataset_record(self, d):
        return {
            "id": d.id, "name": d.name, "description": d.description,
            "createDateTime": self.created, "updateDateTime": self.created, "version": self.data_version,
            "info": {"assemblyId": "GRCh38", "individuals": len(d.biosamples), "genotypedRows": len(d.sites),
                     "variantsObserved": sum(1 for st in d.stats.values() if st.carriers), "synthetic": True},
        }

    def _datasets(self, req, datasets):
        if req.filters or req.params:
            raise BeaconError(400, "filters and request parameters are not supported on datasets in this demo")
        records = [self._dataset_record(d) for d in datasets]
        page = records[req.skip:req.skip + req.limit] if req.limit else records[req.skip:]
        return {"meta": self._meta(req, "record", [spec.returned_schema("dataset")]),
                "responseSummary": {"exists": bool(records), "numTotalResults": len(records)},
                "response": {"collections": page},
                "info": {"notice": NOTICE}}

    # ---- informational endpoints
    def _entry_types(self):
        out = {}
        for entry_type, t in spec.ENTRY_TYPES.items():
            out[entry_type] = {
                "id": entry_type, "name": t["name"], "ontologyTermForThisType": t["ontology"],
                "partOfSpecification": "Beacon v2.0.0", "defaultSchema": spec.schema_reference(entry_type),
                "nonFilteredQueriesAllowed": True,
            }
        out["dataset"]["aCollectionOf"] = [{"id": e, "name": spec.ENTRY_TYPES[e]["name"]}
                                           for e in ("genomicVariant", "individual", "biosample")]
        return out

    def info(self):
        return {"meta": self._info_meta("info", "beacon-info-v2.0.0"),
                "response": {
                    "id": self.beacon_id, "name": self.name, "apiVersion": spec.API_VERSION, "environment": "dev",
                    "organization": self.organization, "version": __version__,
                    "description": ("A small Beacon v2 node serving synthetic genomes: made-up people and genotypes "
                                    "at real GRCh38 positions. For demonstration only; no real people."),
                    "createDateTime": self.created, "updateDateTime": self.created,
                    "info": {"specification": spec.SPEC_TREE, "datasets": list(self.data.datasets)}}}

    def service_info(self):
        return {"id": self.beacon_id, "name": self.name,
                "type": {"group": "org.ga4gh", "artifact": "beacon", "version": spec.API_VERSION},
                "description": "Beacon v2 demo node serving synthetic genomes.",
                "organization": {"name": self.organization["name"],
                                 "url": self.organization.get("welcomeUrl", "https://example.org")},
                "version": __version__, "environment": "dev"}

    def configuration(self):
        return {"meta": self._info_meta("configuration", "beacon-configuration-v2.0.0"),
                "response": {
                    "$schema": f"{spec.SCHEMA_ROOT}/framework/json/configuration/beaconConfigurationSchema.json",
                    "maturityAttributes": {"productionStatus": "DEV"},
                    "securityAttributes": {"defaultGranularity": "boolean",
                                           "securityLevels": list(spec.SECURITY_LEVELS)},
                    "entryTypes": self._entry_types()}}

    def entry_types(self):
        return {"meta": self._info_meta("entryType", "beacon-entry-types-v2.0.0"),
                "response": {"entryTypes": self._entry_types()}}

    def beacon_map(self):
        b = self.base_url
        lookups = {
            "genomicVariant": {"individualLookup": ("individual", "g_variants/{id}/individuals"),
                               "biosampleLookup": ("biosample", "g_variants/{id}/biosamples")},
            "individual": {"genomicVariantLookup": ("genomicVariant", "individuals/{id}/g_variants"),
                           "biosampleLookup": ("biosample", "individuals/{id}/biosamples")},
            "biosample": {"genomicVariantLookup": ("genomicVariant", "biosamples/{id}/g_variants"),
                          "individualLookup": ("individual", "biosamples/{id}/individuals")},
            "dataset": {"genomicVariantLookup": ("genomicVariant", "datasets/{id}/g_variants"),
                        "individualLookup": ("individual", "datasets/{id}/individuals"),
                        "biosampleLookup": ("biosample", "datasets/{id}/biosamples")},
        }
        sets = {}
        for entry_type, related in lookups.items():
            endpoint = spec.ENTRY_TYPES[entry_type]["endpoint"]
            sets[f"{entry_type}Endpoints"] = {
                "entryType": entry_type, "rootUrl": f"{b}/{endpoint}", "singleEntryUrl": f"{b}/{endpoint}/{{id}}",
                "endpoints": {k: {"returnedEntryType": t, "url": f"{b}/{path}"} for k, (t, path) in related.items()},
            }
        return {"meta": self._info_meta("map", "beacon-map-v2.0.0"),
                "response": {"$schema": f"{spec.SCHEMA_ROOT}/framework/json/configuration/beaconMapSchema.json",
                             "endpointSets": sets}}

    def filtering_terms(self):
        terms = [{"type": "ontologyTerm", "id": t, "label": label, "scopes": ["individuals"]}
                 for t, label in sorted(self.data.terms.items())]
        terms += [{"type": "ontologyTerm", "id": t["id"], "label": t["label"], "scopes": ["individuals"]}
                  for t in spec.SEX_TERMS.values()]
        return {"meta": self._info_meta("filteringTerm", "beacon-filtering-terms-v2.0.0"),
                "response": {"filteringTerms": terms, "resources": spec.ONTOLOGY_RESOURCES}}
