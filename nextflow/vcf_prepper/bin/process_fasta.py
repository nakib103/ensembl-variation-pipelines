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
import subprocess
import os
import requests
import glob

from helper import *

FASTA_DIR = "/nfs/production/flicek/ensembl/variation/data/VEP/fasta"

def parse_args(args = None):
    parser = argparse.ArgumentParser()
    
    parser.add_argument(dest="species", type=str, help="species production name")
    parser.add_argument(dest="assembly", type=str, help="assembly default")
    parser.add_argument(dest="version", type=int, help="Ensembl release version")
    parser.add_argument('--division', dest="division", type=str, required = False, help="Ensembl division the species belongs to")
    parser.add_argument('-I', '--ini_file', dest="ini_file", type=str, required = False, help="full path database configuration file, default - DEFAULT.ini in the same directory.")
    parser.add_argument('--fasta_dir', dest="fasta_dir", type=str, required = False, help="FASTA directory")
    parser.add_argument('--force', dest="force", action="store_true")
    
    return parser.parse_args(args)
    
def ungzip_fasta(fasta_dir: str, compressed_fasta: str) -> str:
    if os.path.dirname(compressed_fasta) != fasta_dir:
        print(f"[ERROR] Fasta file {fasta_dir} in wrong directory; should be in - {fasta_dir}")
        exit(1)
        
    process = subprocess.run(["gzip", "-df", compressed_fasta],
        stdout = subprocess.PIPE,
        stderr = subprocess.PIPE
    )
    
    if process.returncode != 0:
        print(f"[ERROR] Could not uncompress fasta file - {compressed_fasta}")
        exit(1)
        
    return compressed_fasta[:-3]
    
def bgzip_fasta(fasta_dir: str, unzipped_fasta: str) -> str:
    if os.path.dirname(unzipped_fasta) != fasta_dir:
        print(f"[ERROR] Fasta file {fasta_dir} in wrong directory; should be in - {fasta_dir}")
        exit(1)
        
    process = subprocess.run(["bgzip", "-f", unzipped_fasta],
        stdout = subprocess.PIPE,
        stderr = subprocess.PIPE
    )
    
    if process.returncode != 0:
        print(f"[ERROR] Could not bgzip fasta file - {unzipped_fasta}")
        exit(1)

    return unzipped_fasta + ".gz"

def index_fasta(zipped_fasta: str, force: str = False) -> None:
    if not os.path.isfile(zipped_fasta):
        print(f"[ERROR] Cannot index fasta - {fasta} - does not exist. Exiting ...")
        exit(1)

    fai = zipped_fasta + ".fai"
    gzi = zipped_fasta + ".gzi"

    if os.path.isfile(fai) and os.path.isfile(gzi) and not force:
        print(f"[INFO] both .fai and .gzi file exist. Skipping ...")
        return

    if os.path.isfile(fai):
        print(f"[INFO] {fai} exist. Deleting ...")
        os.remove(fai)
    
    if os.path.isfile(gzi):
        print(f"[INFO] {gzi} exist. Deleting ...")
        os.remove(gzi)
    
    cmd_index_fasta = "use Bio::DB::HTS::Faidx;"
    cmd_index_fasta += f"Bio::DB::HTS::Faidx->new('{zipped_fasta}');"
    
    process = subprocess.run(["perl", "-e", cmd_index_fasta],
        stdout = subprocess.PIPE,
        stderr = subprocess.PIPE
    )
    if process.returncode != 0:
        print(f"[ERROR] Cannot index fasta file - {zipped_fasta}\n{process.stderr.decode()}\nExiting ...")
        exit(1)
    
def main(args = None):
    args = parse_args(args)
    
    species = args.species
    assembly = args.assembly
    version = args.version
    ini_file = args.ini_file or "DEFAULT.ini"
    core_server = parse_ini(ini_file, "core")
    core_db = get_db_name(core_server, args.version, species, type = "core")
    division = args.division or get_division(core_server, core_db)
    fasta_species_name = get_fasta_species_name(species)

    # TMP - until we use fasta from new website infra
    if species == "homo_sapiens_37":
        fasta_species_name = "Homo_sapiens"
    
    fasta_dir = args.fasta_dir or FASTA_DIR
    fasta_glob = os.path.join(fasta_dir, f"{fasta_species_name}.{assembly}.dna.*.fa.gz")

    fasta = None
    # for human --force won't work; we check and delete fasta manually if needed 
    if glob.glob(fasta_glob) and (not args.force or species.startswith("homo_sapiens")):
        print(f"[INFO] {fasta_glob} exists. Skipping ...")
        
        fasta = os.path.join(fasta_dir, f"{fasta_species_name}.{assembly}.dna.primary_assembly.fa.gz")
        if not os.path.isfile(fasta):
            fasta = os.path.join(fasta_dir, f"{fasta_species_name}.{assembly}.dna.toplevel.fa.gz")
        if not os.path.isfile(fasta):
            print(f"[ERROR] No valid fasta file found, cannot run VEP. Exiting ...")
            exit(1)
    else:
        if glob.glob(fasta_glob):
            print(f"[INFO] {fasta_glob} exists. Will be oerwritten ...")
            for f in glob.glob(fasta_glob):
                os.remove(f)
        
        rl_version = get_relative_version(version, division)
        src_compressed_fasta = get_ftp_path(species, assembly, division, rl_version, "fasta", "local", fasta_species_name)
    
        if src_compressed_fasta is not None:
            compressed_fasta = os.path.join(fasta_dir, os.path.basename(src_compressed_fasta))
            returncode = copyto(src_compressed_fasta, compressed_fasta)
    
        if src_compressed_fasta is None or returncode != 0:
            print(f"[INFO] Failed to copy fasta file - {src_compressed_fasta}, will retry using remote FTP")
        
            compressed_fasta_url = get_ftp_path(species, assembly, division, rl_version, "fasta", "remote", fasta_species_name)
        
            compressed_fasta = os.path.join(fasta_dir, compressed_fasta_url.split('/')[-1])
            returncode = download_file(compressed_fasta, compressed_fasta_url)
            if returncode != 0:
                print(f"[ERROR] Could not download fasta file - {compressed_fasta_url}")
                exit(1)
    
        unzipped_fasta = ungzip_fasta(fasta_dir, compressed_fasta)
        fasta = bgzip_fasta(fasta_dir, unzipped_fasta)
    
    if fasta is not None:
        index_fasta(fasta, force=args.force)
        
if __name__ == "__main__":
    sys.exit(main())
