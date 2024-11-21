#!/usr/bin/env nextflow

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
 
process UPDATE_FIELDS {
  label 'process_medium'

  input: 
  tuple val(meta), path(vcf), path(vcf_index)
  
  output:
  tuple val(meta), path(output_file)
  
  shell:
  output_file =  "UPDATED_S_" + file(vcf).getName()
  source = meta.source
  synonym_file = meta.synonym_file
  rename_clinvar_ids = params.rename_clinvar_ids ? "--rename_clinvar_ids" : ""
  
  '''
  ln -s !{vcf} !{output_file}
  '''
}