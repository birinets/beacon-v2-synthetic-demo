# Changelog

## 0.1.0 (2026-11-18)

First release.

- A GA4GH Beacon v2 node (specification 2.2.0) in the Python standard library, serving a synthetic cohort: 300 made-up people in two datasets, with made-up genotypes at 27 real GRCh38 positions in ten pharmacogenes.
- Sequence, range, bracket and gene queries, with indels matched with or without their anchor base; HPO and sex filters; three access levels (yes or no, counts, records), with counts from 1 to 4 reported as 5; errors instead of answers for unknown parameters, filters, query types and genome builds; and, on every variant answer, whether the queried region was genotyped in each dataset.
- `beacon-demo synth`, `serve`, `query` and `examples`, and a data-directory layout for other data.
- Unit tests, and `tools/validate_schemas.py`, which validates every kind of response and every record in them against the 2.2.0 JSON schemas; on 18 November 2026: 187 responses, 112 records and 6 request bodies, all valid.
