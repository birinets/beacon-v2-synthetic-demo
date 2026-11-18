"""Responses against the Beacon v2.2.0 schemas: required fields, allowed values and types.

The rules below were written by hand from the schema files of https://github.com/ga4gh-beacon/beacon-v2/tree/v2.2.0
(framework/json and models/json/beacon-v2-default-model). They cover what this node returns; they are not a full
JSON Schema validation.
"""

import re
import unittest

from helpers import ask, beacon

CURIE = re.compile(r"^\w[^:]+:.+$")
GRANULARITY = {"boolean", "count", "record"}


def need(test, obj, keys, where):
    test.assertIsInstance(obj, dict, where)
    for key in keys:
        test.assertIn(key, obj, f"{where}: missing {key}")


def ontology_term(test, term, where):
    need(test, term, ["id"], where)
    test.assertRegex(term["id"], CURIE, where)


class Schemas(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.b = beacon()

    def check_request_summary(self, summary):
        need(self, summary, ["apiVersion", "requestedSchemas", "pagination", "requestedGranularity"],
             "receivedRequestSummary")
        self.assertIn(summary["requestedGranularity"], GRANULARITY)
        self.assertIsInstance(summary["requestedSchemas"], list)
        if "filters" in summary:  # Filters: an array of strings
            self.assertTrue(all(isinstance(f, str) for f in summary["filters"]))
        if "includeResultsetResponses" in summary:
            self.assertIn(summary["includeResultsetResponses"], {"ALL", "HIT", "MISS", "NONE"})
        for value in summary.get("requestParameters", {}).values():  # a dictionary of parameter sets
            self.assertIsInstance(value, dict)

    def check_meta(self, meta):  # beaconResponseMeta
        need(self, meta, ["beaconId", "apiVersion", "returnedSchemas", "returnedGranularity", "receivedRequestSummary"],
             "meta")
        self.assertIn(meta["returnedGranularity"], GRANULARITY)
        for schema in meta["returnedSchemas"]:
            need(self, schema, ["entityType", "schema"], "returnedSchemas")
        self.check_request_summary(meta["receivedRequestSummary"])

    def check_entry_response(self, r):
        self.check_meta(r["meta"])
        g = r["meta"]["returnedGranularity"]
        summary = r["responseSummary"]
        need(self, summary, ["exists"], "responseSummary")
        self.assertIsInstance(summary["exists"], bool)
        if g == "count":
            need(self, summary, ["numTotalResults"], "count responseSummary")
            self.assertIn(summary.get("countPrecision", "exact"), {"exact", "imprecise", "rounded"})
            for value in summary.get("countAdjustedTo", [1]):
                self.assertGreaterEqual(value, 1)
        if g == "record":
            need(self, r, ["response"], "resultsets response")
            for result_set in r["response"]["resultSets"]:  # ResultsetInstance
                need(self, result_set, ["id", "setType", "exists"], "result set")
                self.assertIsInstance(result_set["resultsCount"], int)
                for record in result_set.get("results", []):
                    self.assertIsInstance(record, dict)
                    yield_type = r["meta"]["returnedSchemas"][0]["entityType"]
                    getattr(self, f"check_{yield_type}")(record)

    def check_genomicVariant(self, v):
        need(self, v, ["variantInternalId", "variation"], "genomic variation")
        legacy = v["variation"]  # LegacyVariation: variantType, alternateBases and location are required
        need(self, legacy, ["variantType", "alternateBases", "location"], "LegacyVariation")
        for key in ("referenceBases", "alternateBases"):
            self.assertRegex(legacy[key], r"^([ACGTUNRYSWKMBDHV\-\.]*)$")
        location = legacy["location"]  # VRS 1.3 SequenceLocation
        self.assertEqual(location["type"], "SequenceLocation")
        self.assertRegex(location["sequence_id"], CURIE)
        interval = location["interval"]
        self.assertEqual(interval["type"], "SequenceInterval")
        self.assertLess(interval["start"]["value"], interval["end"]["value"])
        for frequency_set in v.get("frequencyInPopulations", []):
            need(self, frequency_set, ["source", "sourceReference", "frequencies"], "FrequencyInPopulations")
            self.assertGreaterEqual(len(frequency_set["frequencies"]), 1)
            for f in frequency_set["frequencies"]:
                need(self, f, ["population", "alleleFrequency"], "PopulationFrequency")
        for case in v.get("caseLevelData", []):
            need(self, case, ["biosampleId"], "CaseLevelVariant")
            ontology_term(self, case["zygosity"], "zygosity")

    def check_individual(self, i):
        need(self, i, ["id", "sex"], "individual")
        ontology_term(self, i["sex"], "sex")
        for feature in i.get("phenotypicFeatures", []):
            ontology_term(self, feature["featureType"], "featureType")

    def check_biosample(self, s):
        need(self, s, ["id", "biosampleStatus", "sampleOriginType"], "biosample")
        ontology_term(self, s["biosampleStatus"], "biosampleStatus")
        ontology_term(self, s["sampleOriginType"], "sampleOriginType")

    def test_entry_responses(self):
        for path, query, level in (
            ("/g_variants", "referenceName=10&start=94781858&referenceBases=G&alternateBases=A", "PUBLIC"),
            ("/g_variants", "geneId=CYP2C9&requestedGranularity=count", "PUBLIC"),
            ("/g_variants", "geneId=CYP2C9&requestedGranularity=record", "CONTROLLED"),
            ("/individuals", "filters=HP:0005110&requestedGranularity=count", "REGISTERED"),
            ("/individuals", "filters=HP:0005110&requestedGranularity=record", "CONTROLLED"),
            ("/biosamples", "requestedGranularity=record&limit=3", "CONTROLLED"),
            ("/g_variants/chr10-94781859-G-A/individuals", "requestedGranularity=record", "CONTROLLED"),
            ("/g_variants/chr16-31096368-C-T/individuals", "filters=HP:0000822&requestedGranularity=count", "REGISTERED"),
        ):
            status, r = ask(self.b, path, query, level)
            self.assertEqual(status, 200, (path, query))
            self.check_entry_response(r)

    def test_collections_response(self):
        status, r = ask(self.b, "/datasets")
        self.check_meta(r["meta"])
        need(self, r, ["responseSummary", "response"], "collections response")
        for dataset in r["response"]["collections"]:
            need(self, dataset, ["id", "name"], "dataset")

    def test_error_response(self):
        status, r = ask(self.b, "/g_variants", "assemblyId=GRCh37&referenceName=1&start=1&referenceBases=A&alternateBases=G")
        self.assertEqual(status, 400)
        self.check_meta(r["meta"])
        need(self, r["error"], ["errorCode"], "BeaconError")
        self.assertIsInstance(r["error"]["errorCode"], int)

    def test_informational_responses(self):
        for path, keys in (("/info", ["id", "name", "apiVersion", "environment", "organization"]),
                           ("/configuration", ["$schema", "maturityAttributes", "entryTypes"]),
                           ("/map", ["$schema", "endpointSets"]),
                           ("/entry_types", ["entryTypes"]),
                           ("/filtering_terms", [])):
            status, r = ask(self.b, path)
            self.assertEqual(status, 200, path)
            need(self, r["meta"], ["beaconId", "apiVersion", "returnedSchemas"], f"{path} meta")
            need(self, r["response"], keys, path)
        status, r = ask(self.b, "/info")
        self.assertIn(r["response"]["environment"], {"prod", "test", "dev", "staging"})
        need(self, r["response"]["organization"], ["id", "name"], "BeaconOrganization")
        status, r = ask(self.b, "/configuration")
        self.assertIn(r["response"]["maturityAttributes"]["productionStatus"], {"DEV", "TEST", "PROD"})
        self.assertTrue(set(r["response"]["securityAttributes"]["securityLevels"]) <= {"PUBLIC", "REGISTERED", "CONTROLLED"})
        for entry_type in r["response"]["entryTypes"].values():
            need(self, entry_type, ["id", "name", "partOfSpecification", "defaultSchema"], "entryTypeDefinition")
            need(self, entry_type["defaultSchema"], ["id", "name", "referenceToSchemaDefinition"], "defaultSchema")
            ontology_term(self, entry_type["ontologyTermForThisType"], "ontologyTermForThisType")
        status, r = ask(self.b, "/map")
        for endpoint_set in r["response"]["endpointSets"].values():
            need(self, endpoint_set, ["entryType", "rootUrl"], "Endpoint")
            for related in endpoint_set.get("endpoints", {}).values():
                need(self, related, ["url", "returnedEntryType"], "RelatedEndpoint")
        status, r = ask(self.b, "/filtering_terms")
        for term in r["response"]["filteringTerms"]:
            need(self, term, ["type", "id"], "FilteringTermInResponse")
        for resource in r["response"]["resources"]:
            need(self, resource, ["id"], "Resource")
        status, r = ask(self.b, "/service-info")  # GA4GH service-info 1.0.0
        need(self, r, ["id", "name", "type", "organization", "version"], "service-info")
        need(self, r["type"], ["group", "artifact", "version"], "service type")
        need(self, r["organization"], ["name", "url"], "service organization")


if __name__ == "__main__":
    unittest.main()
