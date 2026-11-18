# beacon-v2-synthetic-demo

A small [GA4GH Beacon v2](https://docs.genomebeacons.org/) node that answers questions about **synthetic genomes**. It shows how a centre can make a cohort findable without sending its genomes anywhere: the data stays on the node, and the node answers yes or no, a count, or records, depending on who asks.

It follows the Beacon v2 specification, version 2.2.0, and is written with the Python standard library only: no database, no container, no dependencies.

> **Synthetic data, no real people.** The positions and alleles are real (27 GRCh38 positions in ten pharmacogenes, from PharmCAT's public positions file). The 300 people, their sex, phenotypes and genotypes, and the allele frequencies used to draw the genotypes, are all made up from a fixed seed. Phenotypes were drawn independently of genotypes, so the data holds no association of any kind. Frequencies here say nothing about any real population.

## Run it

```bash
git clone https://github.com/birinets/beacon-v2-synthetic-demo
cd beacon-v2-synthetic-demo
PYTHONPATH=src python3 -m beacon_demo serve
```

The node listens on `http://127.0.0.1:5050/api`. Then, in another terminal:

```bash
# Is CYP2C19 rs4244285 (G>A) seen here? Beacon positions are 0-based: VCF position 94781859 is start=94781858.
curl 'http://127.0.0.1:5050/api/g_variants?referenceName=10&start=94781858&referenceBases=G&alternateBases=A'

# How common is it in each dataset? (an allele-frequency answer)
curl 'http://127.0.0.1:5050/api/g_variants?referenceName=10&start=94781858&referenceBases=G&alternateBases=A&requestedGranularity=record'

# How many people with atrial fibrillation (HP:0005110)? The public gets yes or no; a registered user gets a count.
curl -H 'Authorization: Bearer registered-demo' 'http://127.0.0.1:5050/api/individuals?filters=HP:0005110&requestedGranularity=count'
```

Without a server, `PYTHONPATH=src python3 -m beacon_demo query "/g_variants?geneId=CYP2C9&requestedGranularity=count"` prints one answer. `beacon-demo examples out.json` answers a set of example questions and writes them with their responses. `pip install .` adds the `beacon-demo` command.

## What the node does

**Endpoints.** `/info`, `/service-info`, `/configuration`, `/map`, `/entry_types`, `/filtering_terms`, `/datasets`, and the entry types `g_variants`, `individuals` and `biosamples`, each with single-entry lookups (`/g_variants/{id}`) and cross-lookups (`/g_variants/{id}/individuals`, `/individuals/{id}/g_variants`, `/datasets/{id}/individuals` and so on). GET with query parameters and POST with a Beacon request body both work. In a POST body, `requestParameters` may be nested under `g_variant`, as the 2.2.0 framework schema has it, or flat, as many Beacons send them; responses echo them under `g_variant`. `/api` is optional in paths.

**Variant queries.** Sequence queries (`referenceName`, `start`, `referenceBases`, `alternateBases`), range queries (`start` and `end`: any variant that overlaps), bracket queries (two `start` and two `end` values), and gene queries (`geneId`), with `variantType`, `variantMinLength` and `variantMaxLength`. Positions are 0-based, as in Beacon v2. A deletion or insertion matches whether it is written with its anchor base, as a VCF writes it (`start=94949280, GA>G`), or without it (`start=94949281, A>-`): the node compares variants after trimming the bases both alleles share. `alternateBases=N` stands for any single base.

**Filters.** HPO phenotype terms and NCIT sex terms, applied to individuals and combined with AND. `/filtering_terms` lists them.

**Three access levels.** These show the idea; the tokens are public and protect nothing.

| Level | How | People and samples | Variants |
|---|---|---|---|
| Public | no token | yes or no | records with aggregate frequencies only; yes or no when combined with filters on people |
| Registered | `Authorization: Bearer registered-demo` | counts; counts from 1 to 4 are reported as 5 | as public, with counts for filtered questions |
| Controlled | `Authorization: Bearer controlled-demo` | records | records, including which sample carries the variant and its zygosity |

The first two follow the two Beacons of the European Genomic Data Infrastructure (GDI), as its Estonian node describes them: a public allele-frequency Beacon that returns records with frequencies, for Genome of Europe data, and a sensitive Beacon behind a login that returns counts ([GDI Estonia, data flow](https://gdi.ut.ee/en/technical/data-flow/), updated 21 October 2026). That page says the sensitive Beacon is a proof of concept there and that controlled-access requests are not supported. The controlled level here goes further, to show what record-level access would look like.

**It does not answer what it cannot read.** Unknown parameters, unknown filters, other genome builds (`assemblyId` other than GRCh38) and unsupported query types return an error instead of being ignored, because ignoring them would answer a different question. Each variant answer also says, per dataset, whether the queried region was genotyped at all (`info.queriedRegionGenotyped`): `"exists": false` where nothing was genotyped says nothing about anyone. This is an addition in the free-form `info` field, not part of the specification.

## What it does not do

- No real login or permissions. A real node uses an identity service and a data access committee (in the GDI Starter Kit, Life Science Login and REMS) and serves data through htsget with Crypt4GH encryption; none of that is here.
- No ontology reasoning: a filter matches the exact term, and `includeDescendantTerms` is not applied. Similarity other than `exact` is refused.
- No `aminoacidChange`, HGVS (`genomicAlleleShortForm`) or structural-variant queries; no cohorts, runs or analyses.
- Trimming does not left-align: two spellings of an indel in a repeat that differ by a shift do not match.
- Everything is held in memory; it suits thousands of variants, not whole genomes.
- Every kind of response the node gives, and every record in it, validates against the version 2.2.0 JSON schemas with `jsonschema` (`tools/validate_schemas.py`, run in CI). Not checked: the GA4GH VRS 1.3 location object inside variant records (an external schema, stubbed in the check), formats such as date-time and uri, and the official Beacon verifier.

## Your own data

`beacon-demo synth out/` writes the synthetic cohort in the layout the server reads, and `beacon-demo serve --data out/` serves it. The same layout takes other data:

- `datasets.tsv`: dataset id, name, VCF file name, description.
- `samples.tsv`: individual id, biosample id (the VCF column name), dataset id, sex (`female`, `male` or empty), phenotypes (HPO terms separated by `;`).
- `terms.tsv`: id and label of every phenotype term used.
- one VCF per dataset, GRCh38, plain or gzipped, with GT as the first FORMAT field. Multi-allelic rows are split; symbolic alleles are skipped and reported.

Use it only with data you may share in this way. For people's data, that needs consent, ethics approval and a data protection assessment first; this demo is not built for it.

## Tests

```bash
PYTHONPATH=src:tests python3 -m unittest discover -s tests
```

The tests check the query rules against counts read straight from the synthetic VCF, the access levels, the refusals, the HTTP layer, the round trip through files, and the fields the specification requires. The schema check needs the specification and `jsonschema`:

```bash
git clone --depth 1 --branch v2.2.0 https://github.com/ga4gh-beacon/beacon-v2 /tmp/beacon-v2
pip install jsonschema
python3 tools/validate_schemas.py --schemas /tmp/beacon-v2
```

## Contributing

Issues and pull requests are welcome at the source repository, https://github.com/birinets/beacon-v2-synthetic-demo. The copy at https://github.com/indaga-ai/beacon-v2-synthetic-demo is kept identical. See `CHANGELOG.md` for what each version changed.

## Licence and sources

Code and data files: Apache-2.0, see `LICENSE`. Prose documentation (this README): CC BY 4.0, https://creativecommons.org/licenses/by/4.0/.

The positions in `src/beacon_demo/positions.tsv` were selected from `pharmcat_positions.vcf` of [PharmCAT](https://github.com/PharmGKB/PharmCAT) 3.2.0, and each reference allele was checked against GRCh38. Specification: [ga4gh-beacon/beacon-v2](https://github.com/ga4gh-beacon/beacon-v2/tree/v2.2.0). Phenotype terms: [Human Phenotype Ontology](https://hpo.jax.org/). The allele-matching idea follows the reference rule of [genotype-pipeline-checks](https://github.com/birinets/genotype-pipeline-checks).

Written by Ruslan Plakhuta, who runs [Indaga](https://indaga.ai), a product of CryptoSvit, a one-person company in the Netherlands. This repository holds no Indaga code or data.
