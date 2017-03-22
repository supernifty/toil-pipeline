# toil-pipeline
Toil Pipeline

## Dev status
- [18/02/2017] Wrote simple test pipeline based on GATK Best Practise for Germline Sample (with GATK HaplotypeCaller) [src/test-pipeline/pipeline]
- [18/02/2017] Static file (reference file, dbsnp_vcf, etc) and tools executables are defined in a yaml file kept in a folder called 'config' 

## How to run

### Download the assets
```
wget https://swift.rc.nectar.org.au:8888/v1/AUTH_24376b6176a5466b9f20bee02ee1f182/portable-pipelines-0.1.tgz
tar xvfz portable-pipelines-0.1.tgz
```

* The GATK jar is also required.

### Configure the pipeline
```
virtualenv venv
source ./venv/bin/activate
pip install -r requirements.txt
cd src/test_pipeline/example
```

* Edit ../config/toil-test-config.yaml: set the path to wherever the assets were extracted, and to your GATK jar. 
* Edit ./run_toil_test_pipeline.sh: set the paths correctly.

### Run the pipeline

```
./clean.sh
./run_toil_test_pipeline.sh
```

### Ensure the pipeline completed
The test run creates a number of files in Vespa10-2079GL, the final output is Vespa10-2079GL-REDUCED.sorted.markdups.recal.hap.vcf.


## To-do-list
* Share test-data 
* Implement somatic pipeline with GATK Mutect2
* Re-structure pipeline steps for better reusability (currently job chains are declared in steps' functions)
* Remove unused temporary files (use toil Job Temp/Scratch file and store in JobStore instead)
* Add parameters for specifying memory, cpu, etc
* Modify pipeline to support slurm batch submission
