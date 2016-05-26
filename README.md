# mutation_analysis_scripts
day to day activity for mutation data
The scripts are usually the day to day activity that is performed for any exome project for normal/tumor paired samples.
This is usually done after the mutations are called, for normal/tumor pair I prefer using VarScan2 and Mutect.

Am using CNV tools on these data like ADTex , Control-FREEC, VarScan2 CBS, ExomeCNV

For running the rna seq variant pipeline one can simply follow 

`./rna_seq_variant_pipeline.sh  <output_basename> <fastq folder> <output_folder_loc> [cpus]`

Example:

`./rna_seq_variant_pipeline.sh  test_sample /illumina/fastq/sample_folder /path_to/rnaseq/variants/out 12`
