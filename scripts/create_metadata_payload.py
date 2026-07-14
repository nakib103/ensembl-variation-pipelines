#!/usr/bin/env python3

# See the NOTICE file distributed with this work for additional information
# regarding copyright ownership.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

import sys
import configparser
import argparse
import os
import json
import subprocess
import requests
import re
from uuid import UUID
from cyvcf2 import VCF


def parse_args(args=None):
    """Parse command-line arguments for metadata payload creation.

    Args:
        args (list|None): Optional argument list for testing; if None argparse reads from sys.argv.

    Returns:
        argparse.Namespace: Parsed arguments including api_outdir, input_config, endpoint,
            dataset_type and debug flag.
    """
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--pipeline_outdir",
        dest="pipeline_outdir",
        type=str,
        help="path to a vcf prepper pipeline output directory containing api/ and tracks/ subdirectories",
    )
    parser.add_argument(
        "--input_config",
        dest="input_config",
        type=str,
        help="input_config json file used in vcf_prepper",
    )
    parser.add_argument(
        "--endpoint", dest="endpoint", type=str, help="metadata api url"
    )
    parser.add_argument(
        "--dataset_type",
        dest="dataset_type",
        type=str,
        default="all",
        help="dataset type, accepted values: 'short_variants', 'evidence' or 'all'; Default is 'all'",
    )
    parser.add_argument("--debug", dest="debug", action="store_true")

    return parser.parse_args(args)


def is_valid_uuid(uuid: str):
    """Check whether a string is a canonical UUID.

    Args:
        uuid (str): Candidate UUID string.

    Returns:
        bool: True if the string is a valid canonical UUID, False otherwise.
    """
    try:
        uuid_obj = UUID(uuid)
    except ValueError:
        return False
    return str(uuid_obj) == uuid


def get_variant_count(file: str) -> str:
    """Return the number of variant records in a VCF using bcftools.

    Args:
        file (str): Path to the VCF file.

    Returns:
        int|None: Number of records, or None if the count cannot be determined.
    """
    process = subprocess.run(
        ["bcftools", "index", "--nrecords", file],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )

    try:
        return int(process.stdout.decode().strip())
    except Exception as e:
        print(f"""Could not get count from {file}
        {e}""")
        return None


def get_csq_field_index(csq: str, field: str = "Consequence") -> int:
    """Return the index of a specific CSQ field name within the CSQ header string.

    Args:
        csq (str): CSQ header description string (the part after 'Format: ').
        field (str): Field name to locate (default 'Consequence').

    Returns:
        int|None: Index of the named field within the CSQ fields, or None if not found.
    """
    prefix = "Consequence annotations from Ensembl VEP. Format: "
    csq_list = csq[len(prefix) :].split("|")

    for index, value in enumerate(csq_list):
        if value == field:
            return index

    return None


def get_variant_example(file: str, species: str) -> str:
    """Select an example variant identifier from a VCF for the given species.

    Prefers rs699 for human in a specific range, otherwise chooses a missense_variant
    example, and finally falls back to a first-chromosome example.

    Args:
        file (str): Path to the VCF file.
        species (str): Species production name.

    Returns:
        str: A string of the form 'chrom:pos:id' representing an example variant.
    """
    vcf = VCF(file)

    csq_info_description = vcf.get_header_type("CSQ")["Description"].strip('"')
    consequence_idx = get_csq_field_index(csq_info_description, "Consequence")

    # if human, try to find rs699 in 400kbp range
    if species.startswith("homo_sapiens"):
        if '1' in vcf.seqnames:
            for variant in vcf("1:230500000-230900000"):
                if variant.ID == "rs699":
                    chrom = variant.CHROM
                    pos = variant.POS
                    id = variant.ID
                    return f"{chrom}:{pos}:{id}"

    # find a missense_variant
    for variant in vcf:
        csqs = variant.INFO["CSQ"]
        for csq in csqs.split(","):
            consequence = csq.split("|")[consequence_idx]

            if consequence == "missense_variant":
                chrom = variant.CHROM
                pos = variant.POS
                id = variant.ID
                return f"{chrom}:{pos}:{id}"

    # get some random variant if no missense_variant found
    for variant in vcf(f"{vcf.seqnames[0]}:1000"):
        chrom = variant.CHROM
        pos = variant.POS
        id = variant.ID
        return f"{chrom}:{pos}:{id}"


def get_evidence_count(file: str, csq_field: str) -> int:
    """Count variants that have a non-empty value in a specified CSQ subfield.

    Args:
        file (str): Path to the VCF file.
        csq_field (str): CSQ subfield name to test (e.g. 'PHENOTYPES', 'PUBMED').

    Returns:
        int|None: Count of variants with non-empty value in that CSQ field, or None
            if the field is not present in the CSQ header.
    """
    vcf = VCF(file)

    csq_info_description = vcf.get_header_type("CSQ")["Description"].strip('"')
    csq_field_idx = get_csq_field_index(csq_info_description, csq_field)

    if csq_field_idx is None:
        return None

    # find a missense_variant
    count = 0
    for variant in vcf:
        csqs = variant.INFO["CSQ"]
        for csq in csqs.split(","):
            csq_value = csq.split("|")[csq_field_idx]

            if csq_value != "":
                count += 1
                break

    # do not report 0 count
    count = None if count == 0 else count
    return count


def get_source_meta_from_vcf(file: str) -> dict:
    """Extract source metadata from a VCF file's 'source' header.

    Args:
        file (str): Path to the VCF file.

    Returns:
        dict: Parsed key/value pairs from the 'source' header entry.
    """
    vcf = VCF(file)
    source_header = vcf.get_header_type("source")
    vcf.close()

    if not "source" in source_header:
        return {}

    header_content = source_header["source"]
    _, source_info_line = header_content.split('" ', 1)
    source_info = dict(re.findall('(.+?)="(.+?)"\s*', source_info_line))

    return source_info


def parse_input_config(input_config: str) -> dict:
    """Parse an input_config JSON to a mapping keyed by genome UUID.

    Args:
        input_config (str): Path to the input_config JSON file.

    Returns:
        dict|list: Mapping {genome_uuid: {species, assembly}} or an empty list if the file is missing.
    """
    if not os.path.isfile(input_config):
        return []

    with open(input_config, "r") as file:
        input_config_json = json.load(file)

    species_metadata = {}
    for species in input_config_json:
        for genome in input_config_json[species]:
            genome_uuid = genome["genome_uuid"]
            if genome_uuid not in species_metadata:
                species_metadata[genome_uuid] = {}

            species_metadata[genome_uuid]["species"] = genome["species"]
            species_metadata[genome_uuid]["assembly"] = genome["assembly"]
            species_metadata[genome_uuid]["source_name"] = genome["source_name"]
            if "sources" in genome:
                species_metadata[genome_uuid]["sources"] = genome["sources"]

    return species_metadata


def submit_payload(endpoint: str, payload: str) -> str:
    """Submit a payload to a metadata API endpoint via HTTP PUT.

    Args:
        endpoint (str): API endpoint URL.
        payload (str|dict): Payload to send.

    Returns:
        None
    """
    requests.put(endpoint, payload)


def main(args=None):
    """Collect statistics from api VCFs and submit or print payloads.

    Iterates over genomes in the API output directory, computes counts/attributes for
    either 'short_variants' or 'evidence' dataset types, and either prints aggregated payloads
    (debug mode) or submits them to the provided endpoint.

    Args:
        args (list|None): Optional argument list for testing; if None uses sys.argv.

    Returns:
        None
    """
    args = parse_args(args)

    pipeline_outdir = args.pipeline_outdir or os.getcwd()
    input_config = args.input_config or None
    endpoint = args.endpoint or None
    debug = args.debug

    api_outdir = os.path.join(pipeline_outdir, "api")
    tracks_outdir = os.path.join(pipeline_outdir, "tracks")

    if args.dataset_type == "all" or args.dataset_type == None:
        dataset_types = ["short_variants", "evidence"]
    else:
        dataset_types = [args.dataset_type]

    if not debug and endpoint is None:
        print(
            "[ERROR] please provide an endpoint using --endpoint if not using debug mode"
        )
        exit(1)

    species_metadata = {}
    if input_config is not None:
        species_metadata = parse_input_config(input_config)

    for dataset_type in dataset_types:
        if debug:
            aggregate_payload = []

        print(
            f"[INFO] checking directory - {api_outdir} for {dataset_type} statistics data"
        )
        for genome_uuid in os.listdir(api_outdir):
            if species_metadata and genome_uuid not in species_metadata:
                continue

            if not is_valid_uuid(genome_uuid):
                print(f"[WARN] {genome_uuid} is not a valid uuid")
                continue

            # TBD: get this data from thoas if input_config not given
            species = species_metadata[genome_uuid]["species"]
            assembly = species_metadata[genome_uuid]["assembly"]
            source_unparsed = species_metadata[genome_uuid]["source_name"]
            
            source = source_unparsed.replace("%20", " ").replace("%2F", "/")

            api_vcf = os.path.join(api_outdir, genome_uuid, "variation.vcf.gz")
            if not os.path.isfile(api_vcf):
                print(f"[WARN] file not found - {api_vcf}")
                continue

            track_bb = os.path.join(tracks_outdir, genome_uuid, f"variant-{source_unparsed.lower()}-details.bb")
            if not os.path.isfile(track_bb):
                print(f"[WARN] file not found - {track_bb}")
                continue
            track_bw = os.path.join(tracks_outdir, genome_uuid, f"variant-{source_unparsed.lower()}-summary.bw")
            if not os.path.isfile(track_bw):
                print(f"[WARN] file not found - {track_bw}")
                continue
            
            source_meta = get_source_meta_from_vcf(api_vcf)

            payload = {}
            payload["team_responsible"] = "Ensembl Variation"

            payload["name"] = source
            source_version = source_meta.get("version", "")
            if source.lower() == "dbsnp" and not source_version.startswith("build"):
                version = f"build {source_version}"
            elif source.lower() == "eva" and not source_version.startswith("release"):
                version = f"release {source_version}"
            else:
                version = source_version
            payload["label"] = version
            payload["dataset_type"] = dataset_type

            dataset_source = {}
            dataset_source["name"] = api_vcf
            dataset_source["type"] = "vcf"
            payload["dataset_source"] = dataset_source

            source_index = {}
            if os.path.isfile(f"{api_vcf}.tbi"):
                source_index["name"] = f"{api_vcf}.tbi"
            elif os.path.isfile(f"{api_vcf}.csi"):
                source_index["name"] = f"{api_vcf}.csi"
            else:
                print(f"[ERROR] index file not found - {api_vcf}.tbi or {api_vcf}.csi")
                exit(1)
            source_index["type"] = "tabix"
            payload["source_index"] = source_index

            payload["genome_uuid"] = genome_uuid

            dataset_attribute = []

            if dataset_type == "short_variants":
                variant_count = get_variant_count(api_vcf)
                if variant_count is not None:
                    attribute = {}
                    attribute["name"] = "variation.stats.short_variants"
                    attribute["value"] = str(variant_count)
                    dataset_attribute.append(attribute)

                variant_example = get_variant_example(api_vcf, species)
                attribute = {}
                attribute["name"] = "variation.sample_variant"
                attribute["value"] = variant_example
                dataset_attribute.append(attribute)
            else:
                phenotype_count = get_evidence_count(api_vcf, "PHENOTYPES")
                if phenotype_count is not None:
                    attribute = {}
                    attribute["name"] = (
                        "variation.stats.short_variants_with_phenotype_assertions"
                    )
                    attribute["value"] = phenotype_count
                    dataset_attribute.append(attribute)

                publication_count = get_evidence_count(api_vcf, "PUBMED")
                if publication_count is not None:
                    attribute = {}
                    attribute["name"] = (
                        "variation.stats.short_variants_with_publications"
                    )
                    attribute["value"] = publication_count
                    dataset_attribute.append(attribute)

                if species == "homo_sapiens" or species == "homo_sapiens_37":
                    attribute = {}
                    attribute["name"] = (
                        "variation.stats.short_variants_frequency_studies"
                    )
                    attribute["value"] = 1
                    dataset_attribute.append(attribute)

            payload["dataset_attribute"] = dataset_attribute

            payload["track_files"] = {}
            payload["track_files"]["details"] = track_bb
            payload["track_files"]["summary"] = track_bw

            payload["track_source"] = {}
            if source == "MULTIPLE":
                payload["track_source"]["name"] = ", ".join(
                    species_metadata[genome_uuid]["sources"]
                )
            else:
                payload["track_source"]["name"] = source
            payload["track_source"]["url"] = source_meta.get("url", "")
            payload["track_source"]["details"] = version

            if debug:
                aggregate_payload.append(payload)
            else:
                submit_payload(endpoint, payload)

        if debug:
            print(json.dumps(aggregate_payload, indent=4))


if __name__ == "__main__":
    sys.exit(main())
