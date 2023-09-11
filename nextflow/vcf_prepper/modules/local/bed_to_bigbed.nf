#!/usr/bin/env nextflow

/*
* This script generate bigbed file from bed file
*/

process BED_TO_BIGBED {
  input: 
  tuple val(original_vcf), path(bed), val(genome), val(source), val(priorities)
  
  shell:
  output_dir = params.output_dir
  output_bb = file(original_vcf).getName().replace(".vcf.gz", ".bb")
  chrom_sizes = "${projectDir}/assets/chrom_sizes/${genome}.chrom.sizes"
  
  '''
  bedToBigBed -type=bed3+6 !{bed} !{chrom_sizes} !{output_bb}
  
  mkdir -p !{output_dir}/!{genome}/!{source}/tracks
  mv !{output_bb} !{output_dir}/!{genome}/!{source}/tracks/
  '''
}
