"""Validate this node's responses against the GA4GH Beacon v2 JSON schemas, version 2.2.0.

    pip install jsonschema
    git clone --depth 1 --branch v2.2.0 https://github.com/ga4gh-beacon/beacon-v2 /tmp/beacon-v2
    python tools/validate_schemas.py --schemas /tmp/beacon-v2

It sends a fixed set of requests to the node (every example, and each endpoint at each access level and granularity),
validates each response against the framework schema for its kind, every record in it against the default schema
of its entry type, and every request body against the request schema. Exit code 0 when everything validates.

Two substitutions, both stated in the output:

- References to the specification's main branch on raw.githubusercontent.com are read from the local v2.2.0 files.
- The schemas point to GA4GH VRS 1.3 for variant locations, an external schema that is not fetched. Its three
  definitions are stubbed: MolecularVariation and SystemicVariation never match, so a variation must validate as
  the Beacon LegacyVariation, and a VRS Location is accepted as given.

Format keywords such as date-time and uri are not asserted.
"""

from __future__ import annotations

import argparse
import json
import pathlib
import sys
from urllib.parse import urlencode

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "src"))

from jsonschema import Draft202012Validator  # noqa: E402
from referencing import Registry, Resource  # noqa: E402
from referencing.jsonschema import DRAFT202012  # noqa: E402

from beacon_demo import data, synth  # noqa: E402
from beacon_demo.beacon import Beacon  # noqa: E402
from beacon_demo.examples import EXAMPLES, TOKENS  # noqa: E402

MAIN = "https://raw.githubusercontent.com/ga4gh-beacon/beacon-v2/main/"
VRS = "https://w3id.org/ga4gh/schema/vrs/1.3/vrs.json"
VRS_STUB = {"definitions": {"MolecularVariation": False, "SystemicVariation": False, "Location": True}}

FRAMEWORK = "framework/json/"
RESPONSES = {
    "boolean": FRAMEWORK + "responses/beaconBooleanResponse.json",
    "count": FRAMEWORK + "responses/beaconCountResponse.json",
    "record": FRAMEWORK + "responses/beaconResultsetsResponse.json",
    "info": FRAMEWORK + "responses/beaconInfoResponse.json",
    "service-info": FRAMEWORK + "responses/ga4gh-service-info-1-0-0-schema.json",
    "configuration": FRAMEWORK + "responses/beaconConfigurationResponse.json",
    "map": FRAMEWORK + "responses/beaconMapResponse.json",
    "entry_types": FRAMEWORK + "responses/beaconEntryTypesResponse.json",
    "filtering_terms": FRAMEWORK + "responses/beaconFilteringTermsResponse.json",
    "datasets": FRAMEWORK + "responses/beaconCollectionsResponse.json",
    "error": FRAMEWORK + "responses/beaconErrorResponse.json",
    "request": FRAMEWORK + "requests/beaconRequestBody.json",
}
MODEL = "models/json/beacon-v2-default-model/"
RECORDS = {
    "genomicVariant": MODEL + "genomicVariations/defaultSchema.json",
    "individual": MODEL + "individuals/defaultSchema.json",
    "biosample": MODEL + "biosamples/defaultSchema.json",
    "dataset": MODEL + "datasets/defaultSchema.json",
}

VARIANT = "chr10-94781859-G-A"
SEQUENCE = {"referenceName": "10", "start": "94781858", "referenceBases": "G", "alternateBases": "A"}


def registry_for(root):
    resources = []
    for path in sorted(root.glob("framework/json/**/*.json")) + sorted(root.glob("models/json/**/*.json")):
        resource = Resource.from_contents(json.loads(path.read_text(encoding="utf-8")),
                                          default_specification=DRAFT202012)
        relative = path.relative_to(root).as_posix()
        resources += [(path.resolve().as_uri(), resource), (MAIN + relative, resource)]
    resources.append((VRS, Resource.from_contents(VRS_STUB, default_specification=DRAFT202012)))
    return Registry().with_resources(resources)


def requests(beacon):
    """(label, method, path, query, body, level), the examples first."""
    out = [(f"example {x['id']}", x["method"], x["path"], urlencode(x.get("query", {})), x.get("body"), level)
           for x in EXAMPLES for level in x["levels"]]
    person = next(iter(beacon.data.individuals.values()))
    paths = [("/g_variants", SEQUENCE), ("/g_variants", {"geneId": "CYP2C9"}),
             ("/g_variants", {"referenceName": "10", "start": "94761000", "end": "94790000"}),
             (f"/g_variants/{VARIANT}", {}), (f"/g_variants/{VARIANT}/individuals", {}),
             (f"/g_variants/{VARIANT}/biosamples", {}), ("/individuals", {"filters": "HP:0000822"}),
             (f"/individuals/{person.id}", {}), (f"/individuals/{person.id}/g_variants", {}),
             (f"/individuals/{person.id}/biosamples", {}), ("/biosamples", {}),
             (f"/biosamples/{person.biosample_id}", {}), (f"/biosamples/{person.biosample_id}/g_variants", {}),
             (f"/biosamples/{person.biosample_id}/individuals", {}),
             ("/datasets/synthetic-cohort-1/g_variants", {"geneId": "DPYD"}),
             ("/datasets/synthetic-cohort-1/individuals", {}), ("/datasets/synthetic-cohort-2/biosamples", {})]
    for level in TOKENS:
        for granularity in ("boolean", "count", "record"):
            for path, query in paths:
                q = dict(query, requestedGranularity=granularity, limit="3")
                out.append((f"{path} {granularity}", "GET", path, urlencode(q), None, level))
        body = {"meta": {"apiVersion": "v2.2.0"},
                "query": {"requestParameters": {"g_variant": {"referenceName": "10", "start": [94781858],
                                                              "referenceBases": "G", "alternateBases": "A"}},
                          "filters": [{"id": "HP:0000822"}], "requestedGranularity": "record",
                          "includeResultsetResponses": "ALL", "pagination": {"skip": 0, "limit": 2}}}
        out.append(("POST /g_variants with filters", "POST", "/g_variants", "", body, level))
    for path in ("/", "/info", "/service-info", "/configuration", "/map", "/entry_types", "/filtering_terms",
                 "/datasets", "/datasets/synthetic-cohort-2"):
        out.append((f"GET {path}", "GET", path, "", None, "PUBLIC"))
    out += [("unknown endpoint", "GET", "/nothing_here", "", None, "PUBLIC"),
            ("wrong method", "PUT", "/g_variants", "", None, "PUBLIC"),
            ("bad granularity", "GET", "/individuals", "requestedGranularity=everything", None, "PUBLIC")]
    return out


def kind(path, status):
    if status != 200:
        return "error"
    head = [p for p in path.split("/") if p][:1]
    head = head[0] if head else "info"
    if head == "datasets" and path.rstrip("/").count("/") <= 2:
        return "datasets"
    return head if head in RESPONSES and head not in ("datasets",) else "entry"


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--schemas", required=True, type=pathlib.Path,
                        help="a checkout of ga4gh-beacon/beacon-v2 at v2.2.0 (holds framework/json and models/json)")
    args = parser.parse_args(argv)
    root = args.schemas.resolve()
    if not (root / RESPONSES["record"]).is_file():
        parser.error(f"{root} does not hold {RESPONSES['record']}")
    registry = registry_for(root)

    def validator(relative):
        return Draft202012Validator({"$ref": (root / relative).as_uri()}, registry=registry)

    validators = {k: validator(v) for k, v in list(RESPONSES.items()) + list(RECORDS.items())}
    beacon = Beacon(data.from_cohort(synth.make_cohort()))
    failures, counted = [], {"responses": 0, "records": 0, "request bodies": 0}

    def check(name, instance, label):
        errors = sorted(validators[name].iter_errors(instance), key=lambda e: list(e.absolute_path))
        for err in errors[:5]:
            where = "/".join(str(p) for p in err.absolute_path) or "(top)"
            failures.append(f"{label}: {name} at {where}: {err.message[:300]}")
        return not errors

    for label, method, path, query, body, level in requests(beacon):
        token = TOKENS[level]
        headers = {"Authorization": f"Bearer {token}"} if token else {}
        status, response = beacon.handle(method, path, query, body, headers)
        label = f"{label} [{level}, HTTP {status}]"
        if body is not None:
            check("request", body, label + " request body")
            counted["request bodies"] += 1
        what = kind(path, status)
        if what == "entry":
            granularity = response["meta"]["returnedGranularity"]
            check(granularity, response, label)
            entity = response["meta"]["returnedSchemas"][0]["entityType"]
            for result_set in response.get("response", {}).get("resultSets", []):
                for record in result_set.get("results", []):
                    check(entity, record, f"{label} record {result_set['id']}")
                    counted["records"] += 1
        else:
            check(what, response, label)
            for record in response.get("response", {}).get("collections", []) if what == "datasets" else []:
                check("dataset", record, f"{label} dataset")
                counted["records"] += 1
        counted["responses"] += 1

    print(f"Beacon v2 schemas from {root}")
    print("Substitutions: main-branch references read from the local files; VRS 1.3 stubbed "
          "(MolecularVariation and SystemicVariation never match, Location accepted); formats not asserted.")
    print(", ".join(f"{n} {k}" for k, n in counted.items()) + " validated.")
    if failures:
        print(f"{len(failures)} problems:")
        for line in failures:
            print("  " + line)
        return 1
    print("All valid.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
