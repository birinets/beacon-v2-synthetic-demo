import gzip
import json
import os
import shutil
import tempfile
import threading
import unittest
import urllib.error
import urllib.request

from helpers import beacon, cohort

from beacon_demo import data, synth
from beacon_demo.beacon import Beacon
from beacon_demo.server import make_server


class Files(unittest.TestCase):
    def test_the_same_seed_gives_the_same_files(self):
        self.assertEqual(synth.vcf_text(synth.make_cohort(7), "synthetic-cohort-1"),
                         synth.vcf_text(synth.make_cohort(7), "synthetic-cohort-1"))
        self.assertNotEqual(synth.vcf_text(synth.make_cohort(7), "synthetic-cohort-1"),
                            synth.vcf_text(synth.make_cohort(8), "synthetic-cohort-1"))

    def test_a_written_data_directory_loads_to_the_same_index(self):
        tmp = tempfile.mkdtemp()
        try:
            synth.write(cohort(), tmp)
            # a gzipped VCF is read as well
            path = os.path.join(tmp, "synthetic-cohort-2.vcf")
            with open(path, "rb") as src, gzip.open(path + ".gz", "wb") as dst:
                shutil.copyfileobj(src, dst)
            os.remove(path)
            with open(os.path.join(tmp, "datasets.tsv"), encoding="utf-8") as fh:
                rows = fh.read()
            with open(os.path.join(tmp, "datasets.tsv"), "w", encoding="utf-8") as fh:
                fh.write(rows.replace("\tsynthetic-cohort-2.vcf\t", "\tsynthetic-cohort-2.vcf.gz\t"))
            loaded = data.load(tmp)
        finally:
            shutil.rmtree(tmp)
        direct = data.from_cohort(cohort())
        self.assertEqual(loaded.individuals, direct.individuals)
        self.assertEqual(loaded.variants, direct.variants)
        for dataset_id, dataset in direct.datasets.items():
            self.assertEqual(loaded.datasets[dataset_id].stats, dataset.stats)
            self.assertEqual(loaded.datasets[dataset_id].sites, dataset.sites)

    def test_a_vcf_sample_missing_from_the_sample_sheet_stops_loading(self):
        tmp = tempfile.mkdtemp()
        try:
            synth.write(cohort(), tmp)
            sheet = os.path.join(tmp, "samples.tsv")
            with open(sheet, encoding="utf-8") as fh:
                lines = fh.read().splitlines()
            with open(sheet, "w", encoding="utf-8") as fh:
                fh.write("\n".join(line for line in lines if "SMP-0001\t" not in line) + "\n")
            with self.assertRaisesRegex(ValueError, "not in samples.tsv"):
                data.load(tmp)
        finally:
            shutil.rmtree(tmp)


class Server(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.beacon = beacon()
        cls.server = make_server(cls.beacon, port=0)
        cls.url = f"http://127.0.0.1:{cls.server.server_address[1]}"
        cls.thread = threading.Thread(target=cls.server.serve_forever, daemon=True)
        cls.thread.start()

    @classmethod
    def tearDownClass(cls):
        cls.server.shutdown()
        cls.server.server_close()

    def fetch(self, path, body=None, headers=None, method=None):
        data_bytes = body if isinstance(body, bytes) else json.dumps(body).encode() if body is not None else None
        req = urllib.request.Request(self.url + path, data=data_bytes, headers=headers or {}, method=method)
        try:
            with urllib.request.urlopen(req, timeout=10) as resp:
                return resp.status, dict(resp.headers), resp.read()
        except urllib.error.HTTPError as err:
            return err.code, dict(err.headers), err.read()

    def test_http_gives_what_handle_gives(self):
        query = "referenceName=10&start=94781858&referenceBases=G&alternateBases=A&requestedGranularity=record"
        status, headers, body = self.fetch(f"/api/g_variants?{query}")
        self.assertEqual(status, 200)
        self.assertEqual(headers["Access-Control-Allow-Origin"], "*")
        self.assertEqual(json.loads(body), self.beacon.handle("GET", "/g_variants", query)[1])
        status, _, body = self.fetch("/g_variants?" + query)  # without the /api prefix
        self.assertEqual(json.loads(body), self.beacon.handle("GET", "/g_variants", query)[1])

    def test_post_preflight_and_bad_json(self):
        body = {"meta": {"apiVersion": "v2.2.0"},
                "query": {"filters": [{"id": "HP:0005110"}], "requestedGranularity": "count"}}
        status, _, raw = self.fetch("/api/individuals", body, {"Content-Type": "application/json",
                                                                "Authorization": "Bearer registered-demo"})
        self.assertEqual((status, json.loads(raw)["meta"]["returnedGranularity"]), (200, "count"))
        status, headers, _ = self.fetch("/api/individuals", method="OPTIONS")
        self.assertEqual(status, 204)
        self.assertIn("Authorization", headers["Access-Control-Allow-Headers"])
        status, _, raw = self.fetch("/api/individuals", b"{not json", {"Content-Type": "application/json"})
        self.assertEqual(status, 400)
        self.assertEqual(json.loads(raw)["error"]["errorCode"], 400)


class Examples(unittest.TestCase):
    def test_every_example_answers_as_its_text_says(self):
        from beacon_demo.examples import run
        results = {(x["id"], x["accessLevel"]): x for x in run(Beacon(data.from_cohort(cohort())))}
        statuses = {key: x["status"] for key, x in results.items()}
        self.assertEqual({k for k, s in statuses.items() if s != 200},
                         {("wrong-assembly", "PUBLIC"), ("unknown-filter", "PUBLIC")})
        self.assertFalse(results[("genotyped-not-seen", "PUBLIC")]["response"]["responseSummary"]["exists"])
        self.assertFalse(results[("not-genotyped", "PUBLIC")]["response"]["responseSummary"]["exists"])
        self.assertEqual(results[("small-count", "REGISTERED")]["response"]["responseSummary"]["numTotalResults"], 5)
        self.assertEqual(
            results[("deletion-vcf-style", "PUBLIC")]["response"]["responseSummary"],
            results[("deletion-trimmed", "PUBLIC")]["response"]["responseSummary"])


if __name__ == "__main__":
    unittest.main()
