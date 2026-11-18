import unittest

from helpers import ask, beacon, cohort, vcf_counts

RS4244285 = "referenceName=10&start=94781858&referenceBases=G&alternateBases=A"


def record_ids(response):
    return sorted({r["variantInternalId"] for s in response["response"]["resultSets"] for r in s["results"]})


class VariantQueries(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.b = beacon()

    def test_start_is_zero_based(self):
        status, r = ask(self.b, "/g_variants", RS4244285)
        self.assertEqual(status, 200)
        self.assertTrue(r["responseSummary"]["exists"])
        # the VCF position (1-based) by mistake: nothing was genotyped there, and the response says so
        status, r = ask(self.b, "/g_variants", "referenceName=10&start=94781859&referenceBases=G&alternateBases=A")
        self.assertFalse(r["responseSummary"]["exists"])
        self.assertEqual(set(r["info"]["queriedRegionGenotyped"].values()), {False})
        self.assertIn("coverageNote", r["info"])

    def test_genotyped_but_not_seen_is_a_real_no(self):
        status, r = ask(self.b, "/g_variants",
                        "referenceName=13&start=48037781&referenceBases=AGGAGTC&alternateBases=AGGAGTCGGAGTC")
        self.assertFalse(r["responseSummary"]["exists"])
        self.assertEqual(set(r["info"]["queriedRegionGenotyped"].values()), {True})
        self.assertNotIn("coverageNote", r["info"])

    def test_a_deletion_matches_with_or_without_its_anchor_base(self):
        found = []
        for query in ("referenceName=10&start=94949280&referenceBases=GA&alternateBases=G",
                      "referenceName=chr10&start=94949281&referenceBases=A&alternateBases=-",
                      "referenceName=NC_000010.11&start=94949281&referenceBases=A&alternateBases=."):
            status, r = ask(self.b, "/g_variants", query + "&requestedGranularity=record", level="CONTROLLED")
            self.assertEqual(status, 200, r)
            found.append(record_ids(r))
        self.assertEqual(found, [["chr10-94949281-GA-G"]] * 3)

    def test_multi_allelic_rows_are_split_and_n_matches_any_base(self):
        status, r = ask(self.b, "/g_variants", "referenceName=10&start=94942308&referenceBases=G&alternateBases=T"
                                               "&requestedGranularity=record")
        self.assertEqual(record_ids(r), ["chr10-94942309-G-T"])
        status, r = ask(self.b, "/g_variants", "referenceName=10&start=94942308&referenceBases=G&alternateBases=N"
                                               "&requestedGranularity=record")
        self.assertEqual(record_ids(r), ["chr10-94942309-G-A", "chr10-94942309-G-T"])

    def test_range_bracket_and_gene_queries(self):
        status, r = ask(self.b, "/g_variants", "referenceName=10&start=94761000&end=94790000&requestedGranularity=record")
        self.assertEqual(record_ids(r), ["chr10-94761900-C-T", "chr10-94762706-A-G", "chr10-94780653-G-A",
                                         "chr10-94781859-G-A"])
        # a range that ends where a variant starts does not overlap it
        status, r = ask(self.b, "/g_variants", "referenceName=10&start=94781000&end=94781858")
        self.assertFalse(r["responseSummary"]["exists"])
        status, r = ask(self.b, "/g_variants", "referenceName=10&start=94949270,94949290&end=94949280,94949290"
                                               "&requestedGranularity=record")
        self.assertEqual(record_ids(r), ["chr10-94949281-GA-G"])
        status, r = ask(self.b, "/g_variants", "geneId=CYP2C9&requestedGranularity=record&limit=0")
        self.assertEqual(record_ids(r), ["chr10-94942290-C-T", "chr10-94942309-G-A", "chr10-94942309-G-T",
                                         "chr10-94949281-GA-G", "chr10-94981296-A-C", "chr10-94981301-C-G"])
        status, r = ask(self.b, "/g_variants", "geneId=CYP2C9&variantType=DEL&requestedGranularity=record")
        self.assertEqual(record_ids(r), ["chr10-94949281-GA-G"])
        # a gene with a range narrows the range to that gene
        status, r = ask(self.b, "/g_variants", "geneId=CYP2C19&referenceName=10&start=94700000&end=95000000"
                                               "&requestedGranularity=record")
        self.assertEqual(record_ids(r), ["chr10-94761900-C-T", "chr10-94762706-A-G", "chr10-94780653-G-A",
                                         "chr10-94781859-G-A"])

    def test_questions_the_node_cannot_read_are_refused(self):
        for query, words in (
            (RS4244285 + "&assemblyId=GRCh37", "GRCh38 data only"),
            ("aminoacidChange=V600E&geneId=BRAF", "not supported"),
            (RS4244285 + "&datasetIds=x", "unknown query parameter"),
            ("referenceName=10&start=94781858", "needs alternateBases"),
            ("referenceName=25&start=1&referenceBases=A&alternateBases=G", "not a GRCh38 chromosome"),
            ("referenceName=10&start=5&end=5", "greater than start"),
            ("requestedGranularity=everything", "requestedGranularity"),
        ):
            status, r = ask(self.b, "/g_variants", query)
            self.assertEqual(status, 400, query)
            self.assertIn(words, r["error"]["errorMessage"])
        status, r = ask(self.b, "/g_variants", RS4244285 + "&assemblyId=GRCh38.p14")
        self.assertEqual(status, 200)


class Frequencies(unittest.TestCase):
    def test_counts_and_zygosity_agree_with_the_vcf(self):
        b = beacon()
        status, r = ask(b, "/g_variants", RS4244285 + "&requestedGranularity=record", level="CONTROLLED")
        for result_set in r["response"]["resultSets"]:
            ac, an, carriers = vcf_counts(result_set["id"], "chr10", 94781859, "A")
            (record,) = result_set["results"]
            (freq,) = record["frequencyInPopulations"][0]["frequencies"]
            self.assertEqual((freq["alleleCount"], freq["alleleNumber"]), (ac, an))
            self.assertAlmostEqual(freq["alleleFrequency"], ac / an, places=6)
            cases = {c["biosampleId"]: c["zygosity"]["label"] for c in record["caseLevelData"]}
            self.assertEqual(cases, {s: "homozygous" if n == 2 else "heterozygous" for s, n in carriers.items()})
            self.assertEqual(freq["alleleCountHomozygous"], sum(1 for n in carriers.values() if n == 2))

    def test_variant_records_carry_no_people_below_controlled(self):
        b = beacon()
        for level in ("PUBLIC", "REGISTERED"):
            status, r = ask(b, "/g_variants", RS4244285 + "&requestedGranularity=record", level=level)
            self.assertEqual(r["meta"]["returnedGranularity"], "record")
            for result_set in r["response"]["resultSets"]:
                for record in result_set["results"]:
                    self.assertNotIn("caseLevelData", record)


class PeopleAndAccess(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.b = beacon()
        cls.people = cohort().people

    def test_each_level_gets_its_granularity(self):
        expected = sum(1 for p in self.people if "HP:0005110" in p.phenotypes)
        returned = {}
        for level in ("PUBLIC", "REGISTERED", "CONTROLLED"):
            status, r = ask(self.b, "/individuals", "filters=HP:0005110&requestedGranularity=record&limit=0", level)
            returned[level] = r["meta"]["returnedGranularity"]
            if level != "PUBLIC":
                self.assertEqual(r["responseSummary"]["numTotalResults"], expected)
        self.assertEqual(returned, {"PUBLIC": "boolean", "REGISTERED": "count", "CONTROLLED": "record"})
        ids = sorted(x["id"] for s in r["response"]["resultSets"] for x in s["results"])
        self.assertEqual(ids, sorted(p.individual_id for p in self.people if "HP:0005110" in p.phenotypes))

    def test_filters_are_combined_with_and(self):
        status, r = ask(self.b, "/individuals", "filters=HP:0000822,NCIT:C16576&requestedGranularity=count",
                        "CONTROLLED")
        expected = sum(1 for p in self.people if "HP:0000822" in p.phenotypes and p.sex == "female")
        self.assertEqual(r["responseSummary"]["numTotalResults"], expected)

    def test_small_counts_are_rounded_up_below_controlled(self):
        path = "/g_variants/chr16-31096368-C-T/individuals"
        status, exact = ask(self.b, path, "filters=HP:0000822&requestedGranularity=count", "CONTROLLED")
        true_count = exact["responseSummary"]["numTotalResults"]
        self.assertTrue(0 < true_count < 5, "the example needs a count from 1 to 4")
        status, r = ask(self.b, path, "filters=HP:0000822&requestedGranularity=count", "REGISTERED")
        self.assertEqual(r["responseSummary"],
                         {"exists": True, "numTotalResults": 5, "countPrecision": "rounded", "countAdjustedTo": [5]})

    def test_phenotype_filter_on_a_variant_query_hides_the_count_from_the_public(self):
        status, r = ask(self.b, "/g_variants", RS4244285 + "&filters=HP:0000822&requestedGranularity=count")
        self.assertEqual(r["meta"]["returnedGranularity"], "boolean")

    def test_bad_filters_and_tokens_are_refused(self):
        for filters, words in (("HP:0001627", "unknown filter"),):
            status, r = ask(self.b, "/individuals", f"filters={filters}")
            self.assertEqual((status, words in r["error"]["errorMessage"]), (400, True))
        body = {"meta": {"apiVersion": "v2.2.0"}, "query": {"filters": [{"id": "age", "operator": ">", "value": "50"}]}}
        status, r = ask(self.b, "/individuals", method="POST", body=body)
        self.assertEqual(status, 400)
        body["query"]["filters"] = [{"id": "HP:0000822", "scope": "biosamples"}]
        status, r = ask(self.b, "/individuals", method="POST", body=body)
        self.assertEqual(status, 400)
        status, r = self.b.handle("GET", "/individuals", "", None, {"Authorization": "Bearer guess"})
        self.assertEqual(status, 401)
        status, r = ask(self.b, "/no_such_endpoint")
        self.assertEqual(status, 404)
        self.assertIn("meta", r)

    def test_lookups_follow_one_person(self):
        person = next(p for p in self.people if p.dataset_id == "synthetic-cohort-2")
        status, r = ask(self.b, f"/individuals/{person.individual_id}/g_variants", "requestedGranularity=record&limit=0",
                        "CONTROLLED")
        carried = record_ids(r)
        expected = []
        for site, call in zip(cohort().sites, cohort().genotypes[person.biosample_id]):
            for index, alt in enumerate(site.alts, start=1):
                if call and index in call:
                    expected.append(f"{site.chrom}-{site.pos}-{site.ref}-{alt}")
        self.assertEqual(carried, sorted(expected))
        status, r = ask(self.b, f"/biosamples/{person.biosample_id}/individuals", "requestedGranularity=record",
                        "CONTROLLED")
        self.assertEqual([x["id"] for s in r["response"]["resultSets"] for x in s["results"]], [person.individual_id])

    def test_result_sets_follow_include_and_pagination(self):
        query = "referenceName=7&start=99665211&referenceBases=C&alternateBases=T&requestedGranularity=record"
        # rs10264272 was seen in cohort 1 only
        sets = {}
        for include in ("HIT", "MISS", "ALL", "NONE"):
            status, r = ask(self.b, "/g_variants", f"{query}&includeResultsetResponses={include}")
            sets[include] = [s["id"] for s in r["response"]["resultSets"]]
        self.assertEqual(sets, {"HIT": ["synthetic-cohort-1"], "MISS": ["synthetic-cohort-2"],
                                "ALL": ["synthetic-cohort-1", "synthetic-cohort-2"], "NONE": []})
        status, r = ask(self.b, "/individuals", "requestedGranularity=record&skip=2&limit=3", "CONTROLLED")
        first = r["response"]["resultSets"][0]
        self.assertEqual([x["id"] for x in first["results"]], ["IND-0003", "IND-0004", "IND-0005"])
        self.assertEqual(first["resultsCount"], 200)


if __name__ == "__main__":
    unittest.main()
