"""Values taken from the GA4GH Beacon v2 specification, version 2.2.0, and the ontologies it recommends.

The schema identifiers, entry-type names and ontology terms below are copied from the specification's default
model (https://github.com/ga4gh-beacon/beacon-v2/tree/v2.2.0). Nothing here is specific to the synthetic data.
"""

API_VERSION = "v2.2.0"
SPEC_TREE = "https://github.com/ga4gh-beacon/beacon-v2/tree/v2.2.0"
SCHEMA_ROOT = "https://raw.githubusercontent.com/ga4gh-beacon/beacon-v2/v2.2.0"

GRANULARITIES = ("boolean", "count", "record")  # in increasing detail
INCLUDE_RESULTSET_RESPONSES = ("ALL", "HIT", "MISS", "NONE")
SECURITY_LEVELS = ("PUBLIC", "REGISTERED", "CONTROLLED")
DEFAULT_LIMIT = 10

# Entry types this demo serves, as the default model's beaconConfiguration.json and beaconMap.json name them.
ENTRY_TYPES = {
    "genomicVariant": {
        "endpoint": "g_variants",
        "name": "Genomic Variants",
        "ontology": {"id": "ENSGLOSSARY:0000092", "label": "Variant"},
        "schema_id": "ga4gh-beacon-variant-v2.0.0",
        "schema_name": "Default schema for a genomic variation",
        "schema_path": "models/json/beacon-v2-default-model/genomicVariations/defaultSchema.json",
    },
    "individual": {
        "endpoint": "individuals",
        "name": "Individual",
        "ontology": {"id": "NCIT:C25190", "label": "Person"},
        "schema_id": "ga4gh-beacon-individual-v2.0.0",
        "schema_name": "Default schema for an individual",
        "schema_path": "models/json/beacon-v2-default-model/individuals/defaultSchema.json",
    },
    "biosample": {
        "endpoint": "biosamples",
        "name": "Biological Sample",
        "ontology": {"id": "NCIT:C70699", "label": "Biospecimen"},
        "schema_id": "ga4gh-beacon-biosample-v2.0.0",
        "schema_name": "Default schema for a biological sample",
        "schema_path": "models/json/beacon-v2-default-model/biosamples/defaultSchema.json",
    },
    "dataset": {
        "endpoint": "datasets",
        "name": "Dataset",
        "ontology": {"id": "NCIT:C47824", "label": "Data set"},
        "schema_id": "ga4gh-beacon-dataset-v2.0.0",
        "schema_name": "Default schema for datasets",
        "schema_path": "models/json/beacon-v2-default-model/datasets/defaultSchema.json",
    },
}
ENDPOINT_TO_ENTRY_TYPE = {v["endpoint"]: k for k, v in ENTRY_TYPES.items()}

# Recommended terms (individuals and biosamples default schemas; GENO for zygosity).
SEX_TERMS = {
    "female": {"id": "NCIT:C16576", "label": "female"},
    "male": {"id": "NCIT:C20197", "label": "male"},
}
ZYGOSITY_TERMS = {  # (copies of the allele, ploidy of the call) -> GENO term
    (1, 2): {"id": "GENO:0000135", "label": "heterozygous"},
    (2, 2): {"id": "GENO:0000136", "label": "homozygous"},
    (1, 1): {"id": "GENO:0000134", "label": "hemizygous"},
}
BIOSAMPLE_STATUS = {"id": "EFO:0009654", "label": "reference sample"}
SAMPLE_ORIGIN_TYPE = {"id": "OBI:0001479", "label": "specimen from organism"}

ONTOLOGY_RESOURCES = [
    {"id": "hp", "name": "Human Phenotype Ontology", "url": "https://hpo.jax.org/",
     "nameSpacePrefix": "HP", "iriPrefix": "http://purl.obolibrary.org/obo/HP_"},
    {"id": "ncit", "name": "NCI Thesaurus", "url": "https://ncithesaurus.nci.nih.gov/",
     "nameSpacePrefix": "NCIT", "iriPrefix": "http://purl.obolibrary.org/obo/NCIT_"},
]

# GRCh38 RefSeq accessions, for VRS sequence identifiers and HGVS names.
REFSEQ_GRCH38 = {
    "chr1": "NC_000001.11", "chr2": "NC_000002.12", "chr3": "NC_000003.12", "chr4": "NC_000004.12",
    "chr5": "NC_000005.10", "chr6": "NC_000006.12", "chr7": "NC_000007.14", "chr8": "NC_000008.11",
    "chr9": "NC_000009.12", "chr10": "NC_000010.11", "chr11": "NC_000011.10", "chr12": "NC_000012.12",
    "chr13": "NC_000013.11", "chr14": "NC_000014.9", "chr15": "NC_000015.10", "chr16": "NC_000016.10",
    "chr17": "NC_000017.11", "chr18": "NC_000018.10", "chr19": "NC_000019.10", "chr20": "NC_000020.11",
    "chr21": "NC_000021.9", "chr22": "NC_000022.11", "chrX": "NC_000023.11", "chrY": "NC_000024.10",
}


def schema_reference(entry_type):
    t = ENTRY_TYPES[entry_type]
    return {
        "id": t["schema_id"],
        "name": t["schema_name"],
        "referenceToSchemaDefinition": f"{SCHEMA_ROOT}/{t['schema_path']}",
        "schemaVersion": "v2.0.0",
    }


def returned_schema(entry_type):
    return {"entityType": entry_type, "schema": ENTRY_TYPES[entry_type]["schema_id"]}
