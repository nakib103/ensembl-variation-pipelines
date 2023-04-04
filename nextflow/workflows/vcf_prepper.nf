/* Workflow to run Ensembl Variation vcf prepper Pipeline
* The Goal of this workflow is to process and annotate VCF files and generate related track files for genome browser
*/

nextflow.enable.dsl=2

import groovy.json.JsonSlurper
import java.io.File

// params default
params.input_config = "${projectDir}/../nf_config/input_sources.json"
params.output_dir = "/nfs/production/flicek/ensembl/variation/new_website"
params.state = "pre"

// params for nextflow-vep
params.vep_config = "${projectDir}/../nf_config/vep.ini"
params.singularity_dir = "/hps/nobackup/flicek/ensembl/variation/snhossain/website/singularity-images"
params.bin_size = 250000
params.outdir = "/nfs/production/flicek/ensembl/variation/new_website/homo_sapiens_grch38/"

// module imports
repo_dir = "/hps/software/users/ensembl/repositories/${USER}"

include { mergeVCF } from "${projectDir}/../nf_modules/merge_vcf.nf"

include { filterChr } from "${projectDir}/../nf_modules/filter_chr.nf"
include { renameChr } from "${projectDir}/../nf_modules/rename_chr.nf"
include { removeDupIDs } from "${projectDir}/../nf_modules/remove_dup_ids.nf"

include { createFocusVCF } from "${projectDir}/../nf_modules/create_focus_vcf.nf"

include { checkVCF; checkVCF } from "${repo_dir}/ensembl-vep/nextflow/nf_modules/check_VCF.nf"
include { readVCF } from "${repo_dir}/ensembl-vep/nextflow/nf_modules/read_VCF.nf"
include { splitVCF; splitVCF } from "${repo_dir}/ensembl-vep/nextflow/nf_modules/split_VCF.nf"
include { mergeVCF } from "${repo_dir}/ensembl-vep/nextflow/nf_modules/merge_VCF.nf"

include { readChrVCF } from "${projectDir}/../nf_modules/read_chr_vcf.nf"
include { splitChrVCF } from "${projectDir}/../nf_modules/split_chr_vcf.nf"
include { vcfToBed } from "${projectDir}/../nf_modules/vcf_to_bed.nf"
include { bedToBigBed } from "${projectDir}/../nf_modules/bed_to_bigbed.nf"


log.info 'Starting workflow.....'

def slurper = new JsonSlurper()
params.config = slurper.parse(new File(params.input_config))

workflow {
  ch_params = [:]
  ch_params['vcf_file'] = []
  ch_params['index_file'] = []
  ch_params['outdir'] = []
  ch_params['output_prefix'] = []

  genomes = params.config.keySet()
  for (genome in genomes) {
  
    // create a directory for this genome in the output directory
    genome_outdir = params.output_dir + "/" + genome
    file(genome_outdir).mkdir()
  
    for (source in params.config.get(genome)){
      
      // create a directory for this genome in the output directory
      source_outdir = genome_outdir + "/" + source.replace(" ", "_")
      file(source_outdir).mkdir()
      
      // check if index file exist otherwise create it
      index_file = file(source.file_location + ".tbi")
      if (!index_file.exists()){
        def sout = new StringBuilder(), serr = new StringBuilder()
        check_parsing = "tabix -p vcf -f $params.vcf".execute()
        check_parsing.consumeProcessOutput(sout, serr)
        check_parsing.waitFor()
        if( serr ){
          exit 1, "The specified VCF file has issues in parsing: $serr"
        }
      }
      
      // output filename
      prefix = file(source.file_location).getSimpleName() + "_VEP"
      
      ch_params['vcf_file'].add(source.file_location)
      ch_params['index_file'].add(index_file)
      ch_params['outdir'].add(source_outdir)
      ch_params['output_prefix'].add(prefix)
    }
  }
  
  vcf_files = Channel.from(ch_params['vcf_file'])
  index_files = Channel.from(ch_params['index_file'])
  outdir = Channel.from(ch_params['outdir'])

  state = params.state

  //vep(vcf_files, params.vep_config)
  
  //if( state.equals("post")  ) {
    //checkVCF(vcf_files, index_files)
    //readVCF(checkVCF.out, params.bin_size)
    //splitVCF(readVCF.out.transpose())
    //renameChr(splitVCF.out.files.transpose())
    //filterChr(renameChr.out.file)
    //removeDupIDs(filterChr.out.file)
    //mergeVCF(removeDupIDs.out.file.groupTuple())
    //state = "focus"
  //}
  //if ( state.equals("merge") ) {
  //  mergeVCF(ch, source_vcf_outdir)
  //  state = "vep"
  //}
  //if( state.equals("focus")  ) {
  //  createFocusVCF(removeDupIDs.out.vcfFile.collect(), removeDupIDs.out.indexFile.collect(), genome_outdir)
  //  state = "tracks"
  //}
  if( state.equals("tracks")  ) {
    readChrVCF(vcf_files, index_files)
    splitChrVCF(readChrVCF.out.transpose())
    //checkVCF(vcf_files, index_files)
    //readVCF(checkVCF.out, params.bin_size)
    //splitVCF(readVCF.out.transpose())
    vcfToBed(splitChrVCF.out.files.transpose())
    bedToBigBed(vcfToBed.out.bed.groupTuple())
  }
}
