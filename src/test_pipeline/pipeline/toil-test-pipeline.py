'''
Example Toil Pipeline (GATK Best Practise - Germline Pipeline)

Author: Richard Lupat (PeterMac Bioinformatics Core Facility)
Created: 1 Feb 2017
Last Modified: 17 Feb 2017
'''

from argparse import ArgumentParser
import re
import time
import os
import logging
import random
import shutil
import subprocess
import shlex
import sys
import yaml

from toil.job import Job
from toil.common import Toil

class PipelineConfig(object):
	"""
		Store global variables shared by all stages / steps
		e.g. output directory, reference fasta, executables/jars path
	"""
	def __init__(self, outdir, conf_yaml):
		# Check & create output directory
		if (os.path.exists(outdir)):
			print ("This output directory \"%s\" already exists" %(outdir))
			sys.exit(1)
		else:
			try:
				os.makedirs(outdir)
			except:
				print ("Failed to create out directory \"%s\"" %(outdir))
				sys.exit(1)

		self.outdir = os.path.realpath(outdir)

		# Read yaml file, set static file, jars & execs path
		if (os.path.exists(conf_yaml)):
			dat = yaml.load(open(conf_yaml))	
		else:
			print ("This config file \"%s\" does not exist" %(conf_yaml))
			sys.exit(1)
		
		self.static_data = dat['static']
		self.execs = dat['exec']
		self.jars = dat['jar']
			
	def outd(self):
		return self.outdir	

	def get_file(self,k):
		return self.static_data[k]

	def get_exec(self,k):
		return self.execs[k]

	def get_jar(self,k):
		return self.jars[k]

def check_input(job, config, fastqdir):
	'''
		Check Input Fastq Directory
		Find Pair of Fastqs
	'''

	# Input & Output
        ### RL (17/02/2017) Note: Need a more thorough fastq directory & fastq checks
        if os.path.exists(fastqdir):
            fastqs = [ f for f in os.listdir(fastqdir) if os.path.isfile(os.path.join(fastqdir,f)) and f.endswith("fastq.gz")]
        else:
            print("fastq directory '{}' doesn't exist from current directory '{}'".format(fastqdir, os.getcwd()))
            sys.exit(1)

	### RL (17/02/2017)  Temporarily only handle a pair of fastqs
	if len(fastqs) != 2:
		job.fileStore.logToMaster(" Fastq Directory %s does not contain exactly two *fastq.gz files" %(fastqdir), level=logging.CRITICAL)
		return 1
	else:
		if "R1" in fastqs[0]:
			fq1 = os.path.join(fastqdir,fastqs[0])
			fq2 = os.path.join(fastqdir,fastqs[1])
		else:
			fq1 = os.path.join(fastqdir,fastqs[1])
			fq2 = os.path.join(fastqdir,fastqs[0])

		job.fileStore.logToMaster(" Input Fastq1=%s; Fastq2=%s; OutputDirectory=%s " %(fq1, fq2, config.outd()))
		job.addFollowOnJobFn(run_bwa_mem, config, fq1, fq2)

def run_bwa_mem(job, config, fq1, fq2):
	'''
		Run BWA-MEM for alignment against reference Genome
		Input: Pair of FASTQ files
		Output: SAM file
	'''
	# Hardcoded parameters
        ref_fa = config.get_file('ref_fa')

	# BWA
	bwa_exec = config.get_exec('bwa')	
	bwa_mode = "mem"

	# Input & Output 
	samp = re.split(r"[._]+",os.path.basename(fq1))[0]
	outbam = os.path.join(config.outd(),samp + ".sam")

	# Prepare ReadGroup Line
        RG="\"@RG\\tID:%s\\tSM:%s\\tLB:%s\\tPL:ILLUMINA\"" %(samp, samp, samp)

	# Call BWA
	parameters = [bwa_exec, bwa_mode, '-R', RG, '-a', '-M', '-t', '8', ref_fa, fq1, fq2]
	cmd = " ".join(parameters)

	job.fileStore.logToMaster(" Running %s " %(cmd))

	args = shlex.split(cmd)
	fo=open(outbam,"w")
	subprocess.Popen(args, stdout=fo).wait()
	fo.close()

	job.addFollowOnJobFn(sam2bam, config, outbam)

	return outbam

def sam2bam(job, config, inputb):
	'''
		Run SAMTOOLS to convert SAM file to BAM format
		Input: SAM file
		Output: BAM file
	'''

	# Samtools
	samtools_exec = config.get_exec('samtools')
	samtools_mode = "view"
	parameters = "-S -h -b -o"

	# Input & Output
	fn = os.path.splitext(inputb)[0]
	outputbam = fn + ".bam"

	# Call samtools
	cmd = " ".join([samtools_exec, samtools_mode, parameters, outputbam, inputb])
	args = shlex.split(cmd)

	job.fileStore.logToMaster(" Running %s " %(cmd))

	subprocess.Popen(args, stdout=subprocess.PIPE).wait()

	job.addFollowOnJobFn(picard_sortbam, config, outputbam)
	
	return cmd

def picard_sortbam(job, config, inputb):
	'''
		Run PICARD SORTSAM to sort and index bam file
		Input: BAM file
		Output: BAM file (sorted & its bai file)
	'''

        # Helper variables
	tmp_dir = job.fileStore.getLocalTempDir()

        # Hardcoded parameters
        max_record_in_ram = "2000000"   
        create_index = "true"
        validation_stringency = "SILENT"
	sort_order = "coordinate"

        # PICARD MARKDUPS
        picard_jar = config.get_jar('picard')
        picard_exec = "java -Xmx8g -jar %s" %(picard_jar)
        picard_mode = "SortSam"
 
        # Input & Output
        fn = os.path.splitext(inputb)[0]
        outputbam = fn + ".sorted.bam"

        # Call picard markups
        cmd = "%s %s TMP_DIR=%s CREATE_INDEX=%s MAX_RECORDS_IN_RAM=%s SORT_ORDER=%s VALIDATION_STRINGENCY=%s INPUT=%s OUTPUT=%s" %(picard_exec, picard_mode, tmp_dir, create_index, max_record_in_ram, sort_order, validation_stringency, inputb, outputbam)
        args = shlex.split(cmd)

	job.fileStore.logToMaster(" Running %s " %(cmd))

        subprocess.Popen(args, stdout=subprocess.PIPE).wait()

        job.addFollowOnJobFn(picard_markdups, config, outputbam)

        return cmd

def picard_markdups(job, config, inputb):
	'''
		Run PICARD MARKDUPLICATES to mark duplicates reads
		Input: BAM file
		Output: BAM file (sorted & its bai file)
	'''

	# Helper variables
	tmp_dir = job.fileStore.getLocalTempDir()

	# Hardcoded parameters
	max_record_in_ram = "2000000"	
	create_index = "true"
	validation_stringency = "SILENT"

	# PICARD MARKDUPS
	picard_jar = config.get_jar('picard')
	picard_exec = "java -Xmx8g -jar %s" %(picard_jar)
	picard_mode = "MarkDuplicates"

	# Input & Output
	fn = os.path.splitext(inputb)[0]
	metrics = fn + ".metrics.txt"
	outputbam = fn + ".markdups.bam"

	# Call picard markups
	cmd = "%s %s TMP_DIR=%s CREATE_INDEX=%s MAX_RECORDS_IN_RAM=%s METRICS_FILE=%s VALIDATION_STRINGENCY=%s INPUT=%s OUTPUT=%s" %(picard_exec, picard_mode, tmp_dir, create_index, max_record_in_ram, metrics, validation_stringency, inputb, outputbam)
	args = shlex.split(cmd)

	job.fileStore.logToMaster(" Running %s " %(cmd))

	subprocess.Popen(args, stdout=subprocess.PIPE).wait()

	job.addFollowOnJobFn(gatk_base_recal, config, outputbam)

	return cmd

def gatk_base_recal(job, config, inputb):
	'''
		Run GATK BASERECALIBATOR
		Input: Bam file
		Output: .grp file
	'''

	# Hardcoded Parameters
	bed_file	= config.get_file('bed_file')
	ref_fa		= config.get_file('ref_fa')
	known_sites 	= config.get_file('dbsnp_vcf')
	downsampling	= "none"

	# GATK BaseRecalibrator
	gatk_jar = config.get_jar('gatk')
	gatk_exec= "java -Xmx8g -jar %s" %(gatk_jar)
	gatk_mode= "-T BaseRecalibrator"
	
	# Input & Output
	fn = os.path.splitext(inputb)[0]
        output = fn + ".grp"

	# Call GATK BaseRecalibrator
	cmd = "%s %s -I %s -L %s -R %s --downsampling_type %s --knownSites %s -o %s" %(gatk_exec, gatk_mode, inputb, bed_file, ref_fa, downsampling, known_sites, output)
	args = shlex.split(cmd)

	job.fileStore.logToMaster(" Running %s " %(cmd))

        subprocess.Popen(args, stdout=subprocess.PIPE).wait()

	job.addFollowOnJobFn(gatk_printreads, config, output, inputb)

	return cmd

def gatk_printreads(job, config, inputa, inputb):
        '''
                Run GATK PRINTREADS
                Input: .grp file, Bam file
                Output: Bam file (and bai file)
        '''

        # Helper Variables

        # Hardcoded Parameters
        bed_file        = config.get_file('bed_file')
        ref_fa          = config.get_file('ref_fa')        
	downsampling    = "none"

        # GATK PrintReads
        gatk_jar = config.get_jar('gatk')
        gatk_exec= "java -Xmx8g -jar %s" %(gatk_jar)
        gatk_mode= "-T PrintReads"
        
        # Input & Output
        fn = os.path.splitext(inputb)[0]
        output = fn + ".recal.bam"

        # Call GATK PrintReads
	# inputa = .grp
	# inputb = .bam
        cmd = "%s %s --BQSR %s -I %s -L %s -R %s --downsampling_type %s -o %s --filter_bases_not_stored" %(gatk_exec, gatk_mode, inputa, inputb, bed_file, ref_fa, downsampling, output)
        args = shlex.split(cmd)

	job.fileStore.logToMaster(" Running %s " %(cmd))

        subprocess.Popen(args, stdout=subprocess.PIPE).wait()

	job.addFollowOnJobFn(gatk_haplotypecaller, config, output)

        return cmd

def gatk_haplotypecaller(job, config, inputb):
        '''
                Run GATK HAPLOTYPECALLER for variant calling
                Input: Bam file
                Output: VCF file + BAM file
        '''

        # Hardcoded Parameters
        bed_file        	= config.get_file('bed_file')
        ref_fa          	= config.get_file('ref_fa')
	threads			= "1"
	dbsnp			= config.get_file('dbsnp_vcf')
	emitRefConfidence	= "NONE"
        
        # GATK HaplotypeCaller
        gatk_jar = config.get_jar('gatk')
        gatk_exec= "java -Xmx8g -jar %s" %(gatk_jar)
        gatk_mode= "-T HaplotypeCaller"

        # Input & Output
        fn = os.path.splitext(inputb)[0]
        outputbam = fn + ".hap.bam"
	outputvcf = fn + ".hap.vcf"

        # Call GATK HaplotypeCaller
        cmd = "%s %s -I %s -L %s -R %s -o %s --dbsnp %s -nct %s --emitRefConfidence %s -bamout %s" %(gatk_exec, gatk_mode, inputb, bed_file, ref_fa, outputvcf, dbsnp, threads, emitRefConfidence, outputbam)
        args = shlex.split(cmd)

	job.fileStore.logToMaster(" Running %s " %(cmd))

        subprocess.Popen(args, stdout=subprocess.PIPE).wait()

        return cmd

def main():
	parser = ArgumentParser()
	Job.Runner.addToilOptions(parser)

	parser.add_argument('--fastqdir', help='FASTQ Directory', type=str)
	parser.add_argument('--outdir', default='temp', help='OUTPUT Directory', type=str)
	parser.add_argument('--config', default='toil-test-config.yaml', help='Pipeline Config YAML', type=str)
	options = parser.parse_args()

	config = PipelineConfig(options.outdir, options.config)

	print Job.Runner.startToil(Job.wrapJobFn(check_input, config, options.fastqdir),options)

if __name__ == '__main__':
	main()
