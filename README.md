# toil-pipeline
Toil Pipeline

Dev status:
- [18/02/2017] Wrote simple test pipeline based on GATK Best Practise for Germline Sample (with GATK HaplotypeCaller) [src/test-pipeline/pipeline]
- [18/02/2017] Static file (reference file, dbsnp_vcf, etc) and tools executables are defined in a yaml file kept in a folder called 'config' 

How to run:
- Please see [src/test-pipeline/example/run_toil_test_pipeline.sh]

To-do-list:
- Share test-data 
- Implement somatic pipeline with GATK Mutect2
- Re-structure pipeline steps for better reusability (currently job chains are declared in steps' functions)
- Remove unused temporary files (use toil Job Temp/Scratch file and store in JobStore instead)
- Add parameters for specifying memory, cpu, etc
- Modify pipeline to support slurm batch submission
