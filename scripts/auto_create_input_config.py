#!/usr/bin/env python3

"""
Auto-discover EVA-eligible species for the Beta website VCF prepper pipeline.
Generates JSON configs of species requiring re-run based on genebuild status
(or "planned"/"prepared") or EVA variant-count updates.
"""

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

import argparse
import configparser
import gzip
import json
import os
import re
from shutil import rmtree, copy2
import subprocess
import sys
import time
from collections import defaultdict
from pathlib import Path
from typing import List, Optional
import requests

## IMPORTANT: currently this script only supports EVA - we should update the output to have Ensembl data manually -
# - triticum_aestivum
# - triticum_turgidum
# - vitis_vinifera
# - solanum_lycopersicum
# - ovis_aries_rambouillet
# - canis_lupus_familiarisboxer
## These species we can probably ignore until we have EVA data
# - ornithorhynchus_anatinus
# - monodelphis_domestica
# - tetraodon_nigroviridis
## Human will be manual as well


EVA_REST_ENDPOINT = "https://www.ebi.ac.uk/eva/webservices/release"


def parse_args(args=None):
    """
    Parse all command-line arguments for the script.

    Raises:
        SystemExit: Raised internally by *argparse* if `-h/--help`
            is requested or if invalid arguments are supplied.

    Returns:
        argparse.Namespace:  Namespace whose attributes are
            `ini_file`, `output_dir`, and `tmp_dir`, already
            populated with defaults or user-supplied values.
    """

    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument(
        "-I",
        "--ini_file",
        dest="ini_file",
        type=str,
        default="DEFAULT.ini",
        required=False,
        help="Config file with database server information",
    )
    parser.add_argument(
        "-O",
        "--output_dir",
        dest="output_dir",
        type=str,
        default=".",
        required=False,
        help="Dir in which to write output files",
    )

    parser.add_argument(
        "--tmp_dir",
        dest="tmp_dir",
        type=str,
        default=None,
        help="Optional scratch directory for building temporary VCF indexes",
    )

    return parser.parse_args(args)


def parse_ini(ini_file: str, section: str = "database") -> dict:
    """
    Load connection parameters from an INI file.

    Raises:
        SystemExit: If the requested *section* is absent from
            the file (an error message is printed before exit).

    Returns:
        dict:  Keys `host`, `port`, `user` – suitable for passing
            straight to the MySQL command-line client.
    """

    config = configparser.ConfigParser()
    config.read(ini_file)

    if not section in config:
        print(f"[ERROR] Could not find '{section}' config in ini file - {ini_file}")
        exit(1)
    else:
        host = config[section]["host"]
        port = config[section]["port"]
        user = config[section]["user"]

    return {"host": host, "port": port, "user": user}


def get_ensembl_species(server: dict, meta_db: str) -> dict:
    """
    Query metadata db for every available genome.

    Raises:
        SystemExit: If the underlying `mysql` subprocess fails.

    Returns:
        dict:  `{genome_uuid: {species, assembly_id, assembly_name}}`.
    """

    query = f"""
            SELECT
               	g.genome_uuid,
               	g.production_name,
               	a.accession,
               	a.assembly_default,
                o.scientific_name
               	FROM genome AS g, assembly AS a, genome_release AS gr, organism AS o
                WHERE g.assembly_id = a.assembly_id AND g.genome_id = gr.genome_id AND g.organism_id = o.organism_id
                AND suppressed != 1 AND gr.is_current = 1
            """

    process = subprocess.run(
        [
            "mysql",
            "--host",
            server["host"],
            "--port",
            server["port"],
            "--user",
            server["user"],
            "--database",
            meta_db,
            "-N",
            "--execute",
            query,
        ],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )

    if process.returncode != 0:
        print(
            f"[ERROR] Failed to retrieve Ensembl species - {process.stderr.decode().strip()}. \nExiting..."
        )
        exit(1)

    ensembl_species = {}
    for species_meta in process.stdout.decode().strip().split("\n"):
        (genome_uuid, species, assembly, assembly_default, scientific_name) = species_meta.split("\t", 4)
        ensembl_species[genome_uuid] = {
            "species": species,
            "assembly_id": assembly,
            "assembly_name": assembly_default,
            "scientific_name": scientific_name,
        }

    return ensembl_species


def get_ensembl_vcf_filepaths(server: dict, meta_db: str) -> dict:
    """
    Retrieve the location of VCFs already produced by the
    VCF-prepper pipeline.

    Raises:
        SystemExit: On failure of the SQL query / subprocess.

    Returns:
        dict:  `{genome_uuid: {species, assembly_id, file_path}}`.
    """

    query = f"""
            SELECT 
                a.accession,
                g.production_name,
                g.genome_uuid,
                s.name as file_path
            FROM dataset_source AS s
            JOIN dataset AS d 
            ON s.dataset_source_id = d.dataset_source_id
            JOIN dataset_type AS t 
            ON d.dataset_type_id = t.dataset_type_id
            JOIN genome_dataset AS gd 
            ON d.dataset_id = gd.dataset_id
            JOIN genome AS g 
            ON gd.genome_id = g.genome_id
            JOIN assembly AS a
            ON g.assembly_id = a.assembly_id
            WHERE t.name = 'variation'
            AND s.type = 'vcf';
    """

    process = subprocess.run(
        [
            "mysql",
            "--host",
            server["host"],
            "--port",
            server["port"],
            "--user",
            server["user"],
            "--database",
            meta_db,
            "-N",
            "--execute",
            query,
        ],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )

    if process.returncode != 0:
        print(
            f"[ERROR] Failed to retrieve Ensembl species - {process.stderr.decode().strip()}. \nExiting..."
        )
        exit(1)

    ensembl_filepaths = {}
    for filepath_meta in process.stdout.decode().strip().split("\n"):
        (accession, production_name, genome_uuid, file_path) = filepath_meta.split()
        ensembl_filepaths[
            genome_uuid
        ] = {  # using genome_uuid as key, as assembly_id isn't unique i.e. GCA_015227675.2 is mapped to two genome_uuid
            "species": production_name,
            "assembly_id": accession,
            "file_path": file_path,
        }

    return ensembl_filepaths


def get_ensembl_variant_counts(server: dict, meta_db: str) -> dict:
    """
    Pull the pre-computed “short-variant” counts from the
    metadata db.

    Raises:
        SystemExit: When the MySQL subprocess exits non-zero.

    Returns:
        dict:  `{genome_uuid: {"assembly": accession,
                            "variant_count": int}}`.
    """

    query = f"""
            SELECT
                a.accession,
                g.genome_uuid,
                da.value
            FROM dataset_attribute AS da
            JOIN attribute AS attr ON da.attribute_id = attr.attribute_id
            JOIN dataset AS d            ON da.dataset_id       = d.dataset_id
            JOIN dataset_type AS dt      ON d.dataset_type_id   = dt.dataset_type_id
            JOIN genome_dataset AS gd    ON d.dataset_id        = gd.dataset_id
            JOIN genome AS g             ON gd.genome_id        = g.genome_id
            JOIN assembly AS a           ON g.assembly_id       = a.assembly_id
            WHERE
                dt.name = 'variation'
                AND attr.name = 'variation.stats.short_variants';
    """

    process = subprocess.run(
        [
            "mysql",
            "--host",
            server["host"],
            "--port",
            server["port"],
            "--user",
            server["user"],
            "--database",
            meta_db,
            "-N",
            "--execute",
            query,
        ],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )

    if process.returncode != 0:
        print(
            f"[ERROR] Failed to retrieve Ensembl variant counts - {process.stderr.decode().strip()}. \nExiting..."
        )
        exit(1)

    ensembl_variant_counts = {}
    for ensembl_variant_count in process.stdout.decode().strip().split("\n"):
        (assembly, genome_uuid, variant_count) = ensembl_variant_count.split()
        ensembl_variant_counts[genome_uuid] = {
            "assembly": assembly,
            "variant_count": int(variant_count),
        }

    return ensembl_variant_counts


def get_eva_version_from_ensembl_vcf(vcf_path: str):
    """
    Scan the header of the Ensembl VCF and extract
    the EVA release version recorded in the `##source=` line.

    Raises:
        FileNotFoundError: If vcf_path does not exist.
        ValueError:        If the `version="…"` token is present
            but cannot be coerced to an `int`.

    Returns:
        int | None:  The EVA release number, or None when the
            header lacks an EVA `##source` line or version
            isn't an `int`.
    """

    path = Path(vcf_path)
    if not path.exists():
       	new_prefix = "/lts/production/mfreeberg/variation/new_website_production/"
       	new_path = Path("/".join([new_prefix] + vcf_path.split("/")[7:]))
       	path = new_path
       	if not path.exists():
            raise FileNotFoundError(vcf_path)

    opener = gzip.open if path.suffix == ".gz" else open
    src_re = re.compile(r'##source=.*?EVA.*?version="([^"]+)"')

    with opener(path, "rt") as fh:
        for line in fh:
            if not line.startswith("#"):
                break

            if line.startswith("##source=") and "EVA" in line:
                m = src_re.search(line)
                if m:
                    try:
                        return int(m.group(1))
                    except ValueError:
                        return None
                break

    return None


def get_latest_eva_version() -> int:
    """
    Query the EVA REST API for the most recent public release version.

    Raises:
        SystemExit: If the HTTP request fails or an unexpected
            JSON payload is returned.

    Returns:
        int:  The latest EVA release version.
    """

    url = EVA_REST_ENDPOINT + "/v1/info/latest"
    headers = {"Accept": "application/json"}

    response = requests.get(url, headers=headers)

    if response.status_code != 200:
        print(
            f"[ERROR] REST API call for retrieving EVA latest version failed with status - {response.status_code}. Exiting.."
        )
        exit(1)

    try:
        content = response.json()
        release_version = content["releaseVersion"]
    except:
        print(f"[ERROR] Failed to retrieve EVA latest version. Exiting..")
        exit(1)

    return release_version


def get_eva_species(release_version: int) -> dict:
    """
    Fetch per-species statistics for the requested EVA release
    and retain only assemblies with ≥ 5000 RS IDs.

    Raises:
        SystemExit: On HTTP / JSON errors.

    Returns:
        dict:  `{assembly_accession: {species, accession,
                release_folder, taxonomy_id, variant_count}}`.
    """

    eva_species = {}

    url = (
        EVA_REST_ENDPOINT
        + "/v2/stats/per-species?releaseVersion="
        + str(release_version)
    )
    headers = {"Accept": "application/json"}

    response = requests.get(url, headers=headers)
    if response.status_code != 200:
        print(
            f"[ERROR] Could not get EVA species data; REST API call failed with status - {response.status_code}"
        )
        exit(1)

    content = response.json()
    for species in content:
        if species["currentRs"] < 5000:
            continue

        for accession in species["assemblyAccessions"]:
            new_assembly = {
                "species": species["scientificName"],
                "accession": accession,
                "release_folder": species["releaseLink"] or None,
                "taxonomy_id": species["taxonomyId"],
                "variant_count": species["currentRs"],
            }

            eva_species[accession] = eva_species.get(accession, []) + [new_assembly]
    return eva_species


def _tabix_list(path: str, workdir: Path = Path.cwd()) -> Optional[List[str]]:
    """
    Run `tabix -l` and collect the sequence names.

    Raises:
        Any `tabix` error simply yields None.

    Returns:
        list[str] | None:  Sequence names on success; None if
            `tabix` exits non-zero.
    """

    proc = subprocess.run(
        ["tabix", "-l", path],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        cwd=workdir,
    )
    return proc.stdout.strip().splitlines() if proc.returncode == 0 else None


def _header_contigs(path: str) -> List[str]:
    """
    Parse `##contig=<ID=…>` header lines from a VCF
    (compressed or plain).

    Raises:
        IOError:  If the file cannot be read.

    Returns:
        list[str]:  All contig IDs found - in the order they
            appear in the header.
    """

    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt") as fh:
        return [
            ln.split("ID=")[1].split(",")[0].rstrip(">\n")
            for ln in fh
            if ln.startswith("##contig=<ID=")
        ]


def seq_region_matches(eva_file: str, ensembl_file: str, tmp_dir: Path) -> bool:
    """
    Determine whether a remote EVA VCF and a local Ensembl VCF
    reference at least one common sequence region.

    The function attempts `tabix -l` on both files, builds a
    temporary index for the local VCF when necessary, and finally
    falls back to a header parse if indexing fails.

    Raises:
        ValueError:  When *tmp_dir* does not exist.
        Subprocess errors that are thrown while
            copying or indexing files.

    Returns:
        bool:  *True* if a shared contig is found, *False* otherwise.
    """

    # Set up temp_dir/work to hold tmp files
    if not tmp_dir.is_dir():
        raise ValueError(f"--tmp_dir must exist: {tmp_dir}")

    work = tmp_dir / f"work"
    if work.exists():
        print(f"[INFO] Removing stale work directory {work}", file=sys.stderr)
        rmtree(work)
    work.mkdir()

    try:
        # --- EVA (remote) -------------------------------
        max_retries, retry_delay = 4, 60
        eva_seq = None
        for attempt in range(1, max_retries + 1):
            eva_seq = _tabix_list(eva_file, work)
            if eva_seq:
                break
            sys.stderr.write(
                f"[WARN] tabix -l attempt {attempt}/{max_retries} failed "
                f"for {eva_file}\n"
            )
            if attempt < max_retries:
                time.sleep(retry_delay)

        if not eva_seq:
            sys.stderr.write(
                f"[WARN] Cannot obtain seq-regions for remote EVA file "
                f"{eva_file}; skipping comparison.\n"
            )
            return False

        # --- Ensembl local ------------------------------
        ens_seq = _tabix_list(ensembl_file, work)
        if not ens_seq:
            tmp_copy = work / Path(ensembl_file).name
            copy2(ensembl_file, tmp_copy)

            sys.stderr.write(f"[INFO] Building index for {tmp_copy}\n")
            build = subprocess.run(
                ["tabix", "-p", "vcf", "-C", str(tmp_copy)],
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
                cwd=work,
            )

            if build.returncode == 0:
                sys.stderr.write(f"[INFO] Building index successful\n")
                ens_seq = _tabix_list(tmp_copy, work) or []
            else:  # tabix failed – fall back to header parse
                sys.stderr.write(
                    f"[ERROR] Index build failed for {tmp_copy}:\n"
                    f"{build.stderr.strip()}\n"
                )
                try:
                    sys.stderr.write(
                        f"[INFO] Fall-back: attempting header-parse to retrieve contigs\n"
                    )
                    ens_seq = _header_contigs(tmp_copy)
                    sys.stderr.write(
                        f"[INFO] Contigs parsed from header successfully"
                        f"(e.g. {ens_seq[:1]}) ...\n"
                    )
                except Exception as exc:
                    sys.stderr.write(
                        f"[ERROR] Cannot read header of {tmp_copy}: {exc}\n"
                    )
                    return False

        # --- Compare ------------------------------------
        return any(contig in ens_seq for contig in eva_seq)

    # Clean up tmp files
    finally:
        rmtree(work, ignore_errors=True)


def get_ensembl_release_status(server: dict, meta_db: str) -> str:
    """
    Collect “planned” or “prepared” release statuses recorded in
    metadata for every genome.

    Raises:
        SystemExit: If the MySQL query fails.

    Returns:
        dict:  `{genome_uuid: {species, release_status,
                            release_id, assembly_id}}`
    """

    query = f"""
            SELECT 
                g.genome_uuid,
                g.production_name,
                a.accession,
                er.status,
                er.release_id
            FROM genome AS g
            JOIN assembly AS a ON g.assembly_id = a.assembly_id
            JOIN genome_release AS gr ON gr.genome_id    = g.genome_id
            JOIN ensembl_release AS er ON er.release_id  = gr.release_id
            WHERE er.status IN ('prepared','preparing','planned');
    """

    process = subprocess.run(
        [
            "mysql",
            "--host",
            server["host"],
            "--port",
            server["port"],
            "--user",
            server["user"],
            "--database",
            meta_db,
            "-N",
            "--execute",
            query,
        ],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )

    if process.returncode != 0:
        print(
            f"[ERROR] Failed to retrieve release status - {process.stderr.decode().strip()}. \nExiting..."
        )
        exit(1)

    ensembl_release_status = defaultdict(list)
    for release_meta in process.stdout.decode().strip().split("\n"):
        (genome_uuid, production_name, assembly_id, status, release_id) = (
            release_meta.split()
        )
        ensembl_release_status[genome_uuid] = {
            "species": re.sub(
                r"_gca.*$", "", production_name
            ),  # remove accession suffix
            "release_status": status,
            "release_id": release_id,
            "assembly_id": assembly_id,
        }

    return ensembl_release_status


def main(args=None):
    """
    Run the auto-discovery workflow:

    1.  Load CLI arguments and check paths.
    2.  Fetch EVA and Ensembl metadata (species, VCFs, counts, statuses).
    3.  Decide which genomes need to be processed.
    4.  Emit three JSON configuration files.

    Raises:
        SystemExit:  For argument/ environment problems.
        ValueError:  From helpers (e.g. `seq_region_matches`)
            when conditions are not met.

    Returns:
        None
    """

    args = parse_args(args)
    os.makedirs(args.output_dir, exist_ok=True)

    # exit if no tmp dir provided
    if not os.path.isdir(args.tmp_dir):
        sys.stderr.write(
            f"[ERROR] Temporary directory (--tmp_dir) must exist. "
            f"Path given: {args.tmp_dir}\n"
        )
        sys.exit(1)

    # Pull EVA metadata
    eva_release = get_latest_eva_version()
    eva_species = get_eva_species(eva_release)

    # Pull Ensembl metadata
    server = parse_ini(args.ini_file, "metadata")
    ensembl_species = get_ensembl_species(server, meta_db="ensembl_genome_metadata")
    ensembl_vcf_paths = get_ensembl_vcf_filepaths(
        server, meta_db="ensembl_genome_metadata"
    )
    ensembl_status = get_ensembl_release_status(
        server, meta_db="ensembl_genome_metadata"
    )
    ensembl_variant_counts = get_ensembl_variant_counts(
        server, meta_db="ensembl_genome_metadata"
    )

    # Get unique assemblies present in Ensembl
    ensembl_assemblies = [x.get("assembly_id") for x in ensembl_species.values()]
    print(f"[INFO] {len(set(ensembl_assemblies))} unique Ensembl assemblies identified")

    # Get release_id for "Planned" and "Prepared"
    planned_ids = set()
    preparing_ids = set()
    prepared_ids = set()
    for val in ensembl_status.values():
        release_status = val.get("release_status").lower()
        release_id = val.get("release_id")
        if release_status == "planned":
            planned_ids.add(release_id)
        elif release_status == "preparing":
            preparing_ids.add(release_id)
        elif release_status == "prepared":
            prepared_ids.add(release_id)

    # Prepare empty dicts
    ensembl_planned = {}
    ensembl_preparing = {}
    ensembl_prepared = {}
    ensembl_released = {}

    # Loop over only those assemblies present in BOTH Ensembl and EVA
    for asm in set(ensembl_assemblies) & set(eva_species):
        # Get Ensembl genome_uuids associated with assemblies present in BOTH Ensembl and EVA
        associated_uuids = [
            uuid for uuid, rec in ensembl_species.items() if rec["assembly_id"] == asm
        ]

        for uuid in associated_uuids:
            # Grab Ensembl metadata for assembly-uuid
            meta = ensembl_species[uuid]
            status = ensembl_status.get(uuid)
            vcf_meta = ensembl_vcf_paths.get(uuid)

            sp = meta["species"]
            print(
                f"[INFO] Processing... Assembly ID: {asm}. Species: {sp}. Genome: {uuid}."
            )
            scientific_name = meta["scientific_name"]

            # Grab EVA metadata for assembly
            variant_count = -1
            eva_meta = None
            for asm_meta in eva_species[asm]:
                genus = asm_meta["species"].split(" ")[0].lower()
                genus_ensembl = scientific_name.split(" ")[0].lower()
                if genus_ensembl != genus:
                    continue

                if asm_meta["variant_count"] > variant_count:
                    variant_count = asm_meta["variant_count"]
                    eva_meta = asm_meta
            if not eva_meta:
                print(
                    f"[INFO] No EVA metadata found for assembly {asm} with genus {genus_ensembl} matching Ensembl; skipping..."
                )
                continue

            # Skip human
            if sp.startswith("homo"):
                print(f"[INFO] Skipping {sp}")
                continue

            # Build template 'record' dict
            release_folder = eva_meta["release_folder"]
            tax_part = str(eva_meta["taxonomy_id"]) + "_" if eva_release >= 5 else ""
            file_loc = os.path.join(
                release_folder, asm, f"{tax_part}{asm}_current_ids.vcf.gz"
            )

            record = {
                "genome_uuid": uuid,
                "species": meta["species"],
                "assembly": meta["assembly_name"],
                "source_name": "EVA",
                "file_type": "remote",
                "file_location": file_loc,  # EVA file URL
            }

            # Check for genebuild/assembly candidates - candidates must have a planned, preparing, or prepared in the metdata db
            if status:
                print(f"[INFO] Genebuild candidate in the metadata db")
                genome_key = f"{meta['species']}_{meta['assembly_name']}"
                release_id = status.get("release_id", None)
                if release_id is None:
                    print(f"[INFO] Could not find release id for - {uuid}; skipping...")
                    continue

                if status["release_status"].lower() == "prepared":
                    print(f"[INFO] 'prepared' status detected")
                    ensembl_prepared.setdefault(release_id, {}).setdefault(
                        genome_key, []
                    ).append(record)
                if status["release_status"].lower() == "preparing":
                    print(f"[INFO] 'preparing' status detected")
                    ensembl_preparing.setdefault(release_id, {}).setdefault(
                        genome_key, []
                    ).append(record)
                if status["release_status"].lower() == "planned":
                    print(f"[INFO] 'planned' status detected")
                    ensembl_planned.setdefault(release_id, {}).setdefault(
                        genome_key, []
                    ).append(record)

            # Check for EVA variant update candidates that we haven't included it already (as a genebuild/assembly candidate)
            # If we already have a variation data, check if there is any update to EVA
            elif vcf_meta:
                print(
                    f"[INFO] Not genebuild candidate in the metadata db, so checking EVA"
                )
                vcf_path = vcf_meta["file_path"]

                # If the current Ensembl EVA version can be grabbed from vcf header, grab and compare, otherwise revert to variant count comparison
                eva_ensembl_version = get_eva_version_from_ensembl_vcf(
                    vcf_path=vcf_path
                )

                update = False
                if eva_ensembl_version:
                    print(
                        f"[INFO] EVA version parsed from Ensembl file header: {eva_ensembl_version}"
                    )
                    if eva_ensembl_version < eva_release:
                        update = True
                else:
                    print(
                        f"[INFO] EVA version not found in Ensembl file header, comparing variant counts as fall-back"
                    )
                    ensembl_variant_count = ensembl_variant_counts.get(uuid)[
                        "variant_count"
                    ]
                    if ensembl_variant_count < eva_meta["variant_count"]:
                        update = True
                    else:
                        print(
                            f"[INFO] For {sp} EVA variant count has not increased - {ensembl_variant_count} =< {eva_meta['variant_count']}, skipping..."
                        )

                if update:
                    genome_key = f"{meta['species']}_{meta['assembly_name']}"
                    ensembl_released.setdefault(genome_key, []).append(record)

            # If we do not have variation data, it somehow slipped passed previous releases
            else:
                genome_key = f"{meta['species']}_{meta['assembly_name']}"
                ensembl_released.setdefault(genome_key, []).append(record)

    # Write the output JSON files
    print(f"[INFO] Auto-discovery complete. Writing JSON files.")

    for release_id in ensembl_planned:
        planned_json = os.path.join(
            args.output_dir, f"ensembl_planned_{release_id}.json"
        )
        if planned_json:
            with open(planned_json, "w") as fh:
                json.dump(ensembl_planned.get(release_id), fh, indent=4)
        if os.path.isfile(planned_json):
            print(
                f"[INFO] 'Planned' JSON successfully written: {os.path.basename(planned_json)}"
            )

    for release_id in ensembl_preparing:
        preparing_json = os.path.join(
            args.output_dir, f"ensembl_preparing_{release_id}.json"
        )
        if preparing_json:
            with open(preparing_json, "w") as fh:
                json.dump(ensembl_preparing.get(release_id, {}), fh, indent=4)
        if os.path.isfile(preparing_json):
            print(
                f"[INFO] 'Prepared' JSON successfully written: {os.path.basename(preparing_json)}"
            )

    for release_id in ensembl_prepared:
        prepared_json = os.path.join(
            args.output_dir, f"ensembl_prepared_{release_id}.json"
        )
        if prepared_json:
            with open(prepared_json, "w") as fh:
                json.dump(ensembl_prepared.get(release_id, {}), fh, indent=4)
        if os.path.isfile(prepared_json):
            print(
                f"[INFO] 'Prepared' JSON successfully written: {os.path.basename(prepared_json)}"
            )

    released_json = os.path.join(args.output_dir, "ensembl_released.json")
    if ensembl_released:
        with open(released_json, "w") as fh:
            json.dump(ensembl_released, fh, indent=4)
    if os.path.isfile(released_json):
        print(
            f"[INFO] 'Released' JSON successfully written: {os.path.basename(released_json)}"
        )

    print(f"Done.")


if __name__ == "__main__":
    sys.exit(main())
