/*
 * See the NOTICE file distributed with this work for additional information
 * regarding copyright ownership.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
 
//
// Workflow to run Ensembl Variation vcf prepper Pipeline
// The Goal of this workflow is to process and annotate VCF files and generate related track files for genome browser
//

import groovy.json.JsonSlurper
import java.io.File

def slurper = new JsonSlurper()
params.config = slurper.parse(new File(params.input_config))

include { CREATE_RANK_FILE } from "../modules/local/create_rank_file.nf"
include { PREPARE_GENOME } from "../subworkflows/local/prepare_genome.nf"
include { UPDATE_FIELDS } from "../modules/local/update_fields.nf"
include { REMOVE_VARIANTS } from "../modules/local/remove_variants.nf"
include { RUN_VEP } from "../subworkflows/local/run_vep.nf"
include { COUNT_VCF_VARIANT } from "../modules/local/count_vcf_variant.nf"
include { SPLIT_VCF } from "../subworkflows/local/split_vcf.nf"
include { VCF_TO_BED } from "../modules/local/vcf_to_bed.nf"
include { CONCAT_BEDS } from "../modules/local/concat_beds.nf"
include { BED_TO_BIGBED } from "../modules/local/bed_to_bigbed.nf"
include { BED_TO_WIG } from "../modules/local/bed_to_wig.nf"
include { WIG_TO_BIGWIG } from "../modules/local/wig_to_bigwig.nf"
include { SUMMARY_STATS } from "../modules/local/summary_stats.nf"

def parse_config (config) {
  input_set = []
  
  genomes = config.keySet()
  for (genome in genomes) {
    for (source_data in params.config.get(genome)) {
      vcf = source_data.file_location
      
      meta = [:]
      meta.genome = genome
      meta.genome_uuid = source_data.genome_uuid
      meta.species = source_data.species
      meta.assembly = source_data.assembly
      meta.source = source_data.source_name.replace(" ", "_")
      meta.file_type = source_data.file_type
      
      input_set.add([meta, vcf])
    }  
  }
  
  return input_set
}

workflow VCF_PREPPER {
  if (params.skip_vep && params.skip_tracks && params.skip_stats) {
    exit 0, "Skipping VEP and track file generation, nothing to do ..."
  }
  
  input_set = parse_config(params.config)
  ch_input = Channel.fromList( input_set )
  
  // setup
  CREATE_RANK_FILE( params.rank_file )
  PREPARE_GENOME( ch_input )
    
  // api files
  if (!params.skip_vep) {
    // pre-process
    UPDATE_FIELDS( PREPARE_GENOME.out )
    REMOVE_VARIANTS( UPDATE_FIELDS.out )
    
    // run vep
    vep = RUN_VEP( REMOVE_VARIANTS.out )

    // post-process
    COUNT_VCF_VARIANT( vep )

    COUNT_VCF_VARIANT.out
    .map {
      meta, vcf, vcf_index, variant_count ->
        // if vcf has no variant - remove output directories and filter channel 
        if ( variant_count.equals("0") ) {
          file(meta.genome_api_outdir).delete()
          file(meta.genome_tracks_outdir).delete()

          "NO_VARIANT"
        }
        else {
          [meta, vcf, vcf_index]
        }
    }
    .filter { ! it.equals("NO_VARIANT") }
    .set { ch_post_api }
  }
  else {
    ch_post_api = PREPARE_GENOME.out
  }

  // track files
  if (!params.skip_tracks) {
    // create bed from VCF
    // TODO: vcf_to_bed maybe faster without SPLIT_VCF - needs benchmarking
    SPLIT_VCF( ch_post_api )
    VCF_TO_BED( CREATE_RANK_FILE.out, SPLIT_VCF.out.transpose() )
    CONCAT_BEDS( VCF_TO_BED.out.groupTuple() )

    // create source tracks
    // TODO: remove symlink creation for focus track when we have multiple source
    BED_TO_BIGBED( CONCAT_BEDS.out )
    BED_TO_WIG( CONCAT_BEDS.out )
    WIG_TO_BIGWIG( BED_TO_WIG.out )

    // if track generation is run vep-ed VCF file move needs to wait for this step to finish
    SPLIT_VCF.out
    .map {
      meta, splits ->
        [meta]
    }
    .set { ch_split_finish }
  }

  // summary stats
  if (!params.skip_stats) {
    // it can run in parallel to track generation
    SUMMARY_STATS( ch_post_api )

    SUMMARY_STATS.out
    .set { ch_stats_finish }
  }

  // post process
  if (!params.skip_vep || !params.skip_stats){
    if(!params.skip_stats && !params.skip_tracks) {
      ch_split_finish
      .join ( ch_stats_finish )
      .map {
        meta, vcf, vcf_index ->
          [meta, vcf, vcf_index]
      }
      .set { ch_post_process }
    }
    else if (params.skip_stats && !params.skip_tracks) {
      ch_split_finish
      .join ( ch_post_api )
      .map {
        meta, vcf, vcf_index ->
          [meta, vcf, vcf_index]
      }
      .set { ch_post_process }
    }
    else if (!params.skip_stats && params.skip_tracks) {
      ch_stats_finish
      .set { ch_post_process }
    }
    else {
      ch_post_api
      .set { ch_post_process }
    }

    ch_post_process
    .map {
      meta, vcf, vcf_index ->
        // TODO: when we have multiple source per genome we need to delete source specific files
        new_vcf = "${meta.genome_api_outdir}/${vcf}"
        new_vcf_index = "${meta.genome_api_outdir}/${vcf_index}"
        
        // in -resume vcf and vcf_index may not exists as already renamed
        // moveTo instead of renameTo - in -resume dest file may exists from previous run
        if ( file(vcf).exists() && file(vcf_index).exists() ) {
          file(vcf).moveTo(new_vcf)
          file(vcf_index).moveTo(new_vcf_index)
        }

        [meta, new_vcf, new_vcf_index]
    }
  }
}
