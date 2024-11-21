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
from cyvcf2 import VCF, Writer
from Bio import bgzf
import argparse

def parse_args(args = None, description: bool = None):
    parser = argparse.ArgumentParser(description = description)
    
    parser.add_argument(dest="input_file", type=str, help="input VCF file")
    parser.add_argument(dest="source", type=str, help="input VCF file source")
    parser.add_argument('--rename_clinvar_ids', dest="rename_clinvar_ids", action="store_true")
    parser.add_argument('--chromosomes', dest="chromosomes", type=str, help="comma separated list of chromosomes to put in header")
    parser.add_argument('-O', '--output_file', dest="output_file", type=str)
    
    return parser.parse_args(args)

def format_clinvar_id(id: str) -> str:
    if not id.startswith("VCV"):
        leading_zero = 9 - len(id)
        return "VCV" + ("0" * leading_zero) + id

    return id

def main(args = None):
    args = parse_args(args)

    input_file = args.input_file
    source = args.source
    chromosomes = args.chromosomes or None
    output_file = args.output_file or os.path.join(os.path.dirname(input_file), "UPDATED_S_" + os.path.basename(input_file))

    if args.rename_clinvar_ids and source == "ClinVar":
        format_id = format_clinvar_id
    else:
        format_id = lambda x : x

    input_vcf = VCF(input_file)
    output_vcf = Writer(output_file, input_vcf, mode="wz")
    for variant in input_vcf:
        variant.ID = format_id(variant.ID)
        try:
            del variant.INFO["vep"]
        except:
            pass

        output_vcf.write_record(variant)
    input_vcf.close()
    
if __name__ == "__main__":
    sys.exit(main())
