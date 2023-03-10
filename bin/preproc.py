#!/usr/bin/python3
 
"""
Copyright 2023 O. Sotolongo <asqwerty@gmail.com>
 
This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.
 
This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
"""

import sys
import os
import getopt
import re
import subprocess
import gzip
from slurm import send_sbatch

"""
See: 
 - For WES pipeline: http://detritus.fundacioace.com/wiki/doku.php?id=genetica:wes
 - For execution into SLURM: https://github.com/asqwerty666/slurm-modpy
"""

"""
Data PATHS
"""
src_dir = '/ruby/Traspaso-datos-11-2022/fqdata'
ref_dir = '/ruby/WES_reference/hg38'
ref_name = 'Homo_sapiens_assembly38'
ref_fa = ref_dir+'/'+ref_name+'.fasta'
output_dir = '/ruby/WES_output'
tmp_dir = '/ruby/'+os.environ.get('USER')+'/tmp/'
baits = ref_dir + '/KAPA_HyperExome_hg38_bait.interval_list'
targets = ref_dir + '/KAPA_HyperExome_hg38_target.interval_list'
unions = ref_dir + '/KAPA_HyperExome_hg38_union.interval_list' 
known1 = 'Homo_sapiens_assembly38.known_indels.vcf.gz'
known2 = 'Mills_and_1000G_gold_standard.indels.hg38.vcf.gz'
dbsnp = 'Homo_sapiens_assembly38.dbsnp138.vcf'
hcsnps = '1000G_phase1.snps.high_confidence.hg38.vcf.gz'
"""
Executable PATHS
"""
fastqc = '/nas/usr/local/bin/fastqc'
bwa = '/nas/usr/local/bin/bwa mem -t 4 -M'
picard = 'java -Djava.io.tmpdir='+tmp_dir+' -Xmx8g -jar /nas/usr/local/bin/picard.jar'
samtools = '/nas/software/samtools/bin/samtools'
verifyBamID = '/nas/usr/local/bin/verifyBamID'
gatk = 'singularity run --cleanenv -B /nas:/nas -B /ruby:/ruby -B /the_dysk:/the_dysk /nas/usr/local/opt/gatk4.simg gatk --java-options "-DGATK_STACKTRACE_ON_USER_EXCEPTION=true -Xmx16G"'
snpEff = 'java -Xmx8g -jar /nas/software/snpEff/snpEff.jar'
deepvariant = 'singularity run --cleanenv -B /nas:/nas -B /ruby:/ruby -B /the_dysk:/the_dysk /nas/usr/local/opt/deepvariant.simg /opt/deepvariant/bin/run_deepvariant'
# Get CLI inputs
short_args = 'c:o:gs:'
long_args = ['cut=', 'output=', 'debug', 'source=']
debug = 0
cfile=''
outdatadir=''
try:
	args, values = getopt.getopt(sys.argv[1:], short_args, long_args)
except getopt.error as err:
	print(str(err))
	sys.exit(2)

for a,v in args:
	if a in ('--cut', '-c'):
		cfile = v
	elif a in ('--output', '-o'):
		outdatadir = v
	elif a in ('--debug', '-g'):
		debug = 1
	elif a in ('--source', '-s'):
		src_dir = v

lpath = os.path.abspath(src_dir)
dir_cont = next(os.walk(src_dir))[1]

if cfile and os.path.isfile(cfile):
	f = open(cfile, 'r')
	cuts = f.read().splitlines()
	dir_cont = list(filter(lambda x: x in set(dir_cont), set(cuts)))

# Creo el entorno de ejecucion
tmp_shit = os.environ.get('TMPDIR')
if not os.path.isdir(tmp_shit):
	tmp_shit = os.environ.get('PWD') + '/tmp/'
	os.mkdir(tmp_shit)
if outdatadir:
	if not os.path.isdir(outdatadir): os.mkdir(outdatadir)
	wdir = outdatadir
else:
	wdir = os.environ.get('PWD')
outdir = wdir+'/slurm'
if not os.path.isdir(outdir): os.mkdir(outdir)
fwdir = wdir+'/final'
if not os.path.isdir(fwdir): os.mkdir(fwdir)
fq = {}
jobids = []
cdata = {'time':'24:0:0', 'cpus':'8', 'mem-per-cpu':'4G'}
pollos_list = wdir+'/subjects.list'
if os.path.exists(pollos_list):  os.remove(pollos_list)
for pollo in dir_cont:
	fq_list = []
	fq_id =''
	fq_files = next(os.walk(lpath+'/'+pollo))[2]
	udir = wdir+'/'+pollo
	if not os.path.isdir(udir): os.mkdir(udir)
	tmpdir = udir+'/tmp'
	if not os.path.isdir(tmpdir): os.mkdir(tmpdir)
	resdir = udir+'/results'
	if not os.path.isdir(resdir): os.mkdir(resdir)
	fjobids = []
	ljobids = []
	gsconv = ''
	with open(pollos_list, 'a') as lf:
		lf.write(pollo+'\t'+resdir + '/' + pollo +'_raw.snps.indels.g.vcf.gz\n')
	# First job, fastqc
	cdata['job_name'] = pollo + '_fastqc'
	cdata['filename'] = outdir + '/' + pollo + '_fastqc.sh'
	cdata['output'] = outdir + '/' + pollo + '_fastqc.out'
	cdata['command'] = 'mkdir -p ' + udir + '/qc\n'
	cdata['command'] += fastqc + ' -o ' + udir + '/qc ' + lpath + '/' + pollo + '/*.fq.gz'
	cdata.pop('dependency', None)
	p = send_sbatch(cdata)
	ljobids.append(p)
	
	for fq in range(2): #Esta es la logica que hay que cambiar cada vez que nos cambian la estructura de los datos
		gsconv += ' -I ' + tmpdir + '/' + pollo + '_' + str(fq+1) + '.sam'
		# Now, bwa for each fq.gz
		cdata['filename'] = outdir + '/' + pollo + '_bwa_' + str(fq) + '.sh'
		cdata['job_name'] = pollo + '_bwa_' + str(fq)
		cdata['output'] = outdir + '/' + pollo + '_bwa_' + str(fq) + '.out'
		cdata['command'] = bwa + ' -R "@RG\\tID:1\\tDS:KAPA_TE\\tPL:ILLUMINA\\tLB:'+pollo+'\\tSM:'+pollo+'" ' + ref_fa+' ' + lpath + '/' + pollo + '/' + pollo + '_'+ str(fq+1) + '.fq.gz > '+ tmpdir + '/' + pollo + '_' + str(fq+1) + '.sam'
		p = send_sbatch(cdata)
		fjobids.append(p)
	
	# Next batch, MergeSamFiles +  SortSam + index
	deps = ',afterok:'.join(map(str, fjobids))
	fjobids = []
	cdata['job_name'] = pollo + '_MergeSortIndex'
	cdata['filename'] = outdir + '/' + pollo + '_MergeSortIndex.sh'
	cdata['output'] = outdir + '/' + pollo + '_MergeSortIndex.out'
	cdata['command'] = gatk + ' MergeSamFiles' + gsconv + ' -O ' + tmpdir + '/' + pollo + '.sam\n'
	cdata['command'] += gatk + ' SortSam -I ' + tmpdir + '/' + pollo + '.sam -O ' + tmpdir + '/' + pollo + '_sorted.bam --SORT_ORDER coordinate\n'
	cdata['command'] += samtools + ' index ' + tmpdir + '/' + pollo + '_sorted.bam'
	cdata['dependency'] = 'afterok:' + deps
	p = send_sbatch(cdata)
	#print(cdata)
	# VerifyBamID for the sorted bam file
	cdata['job_name'] = pollo + '_verifybamid'
	cdata['filename'] = outdir + '/' + pollo + '_verifybamid.sh'
	cdata['output'] = outdir + '/' + pollo + '_verifybamid.out'
	cdata['command'] = verifyBamID + ' --vcf ' + ref_dir + '/hapmap_3.3.hg38.vcf.gz --bam ' + tmpdir + '/' + pollo + '_sorted.bam --chip-none --maxDepth 1000 --precise --verbose --ignoreRG --out ' + resdir + '/' + pollo + '_verifybam |& grep -v "Skipping marker"'
	cdata['dependency'] = 'afterok:' + str(p)
	p0 = send_sbatch(cdata)
	ljobids.append(p0)
	# ValidateSamFile
	cdata['job_name'] = pollo + '_Validate'
	cdata['filename'] = outdir + '/' + pollo + '_Validate.sh'
	cdata['output'] = outdir + '/' + pollo + '_Validate.out'
	cdata['command'] = gatk + ' ValidateSamFile --IGNORE MATE_NOT_FOUND --IGNORE MISSING_READ_GROUP --IGNORE RECORD_MISSING_READ_GROUP -I ' + tmpdir + '/' + pollo + '_sorted.bam'
	cdata['dependency'] = 'afterok:' + str(p)
	p0 = send_sbatch(cdata)
	ljobids.append(p0)
	# MarkDuplicates
	cdata['job_name'] = pollo + '_MarkDuplicates'
	cdata['filename'] = outdir + '/' + pollo + '_MarkDuplicates.sh'
	cdata['output'] = outdir + '/' + pollo + '_MarkDuplicates.out'
	cdata['command'] = gatk +  ' MarkDuplicates -I ' + tmpdir + '/' + pollo + '_sorted.bam -O ' + resdir + '/' + pollo + '_rmdups.bam --METRICS_FILE ' + resdir + '/' + pollo + '_metrics.txt --QUIET TRUE --MAX_RECORDS_IN_RAM 2000000 --ASSUME_SORTED TRUE --CREATE_INDEX TRUE'
	cdata['dependency'] = 'afterok:' + str(p)
	p = send_sbatch(cdata)
	# CollectHsMetrics para sacar un report
	cdata['job_name'] = pollo + '_CollectHsMetrics_1'
	cdata['filename'] = outdir + '/' + pollo + '_CollectHsMetrics_1.sh'
	cdata['output'] = outdir + '/' + pollo + '_CollectHsMetrics_1.out'
	cdata['command'] = gatk + ' CollectHsMetrics -BI ' + baits + ' -TI ' + targets + ' -I ' + resdir + '/' + pollo + '_rmdups.bam -O ' + resdir + '/' + pollo + '_hsmetrics.txt'
	cdata['dependency'] = 'afterok:' + str(p)
	p0 = send_sbatch(cdata)
	ljobids.append(p0)
	# Vamos a GATK4, BaseRecalibrator, ApplyBQSR, AnalyzeCovariates, BaseRecalibrator, AnalyzeCovariates, HaplotypeCaller
	# BaseRecalibrator
	cdata['job_name'] = pollo + '_BaseRecalibrator_1'
	cdata['filename'] = outdir + '/' + pollo + '_BaseRecalibrator_1.sh'
	cdata['output'] = outdir + '/' + pollo + '_BaseRecalibrator_1.out'
	cdata['command'] = gatk + ' BaseRecalibrator -I ' + resdir + '/' + pollo + '_rmdups.bam -R ' + ref_fa +' --known-sites ' + ref_dir + '/'+known1+' --known-sites '+ ref_dir + '/'+known2+' --known-sites ' + ref_dir + '/'+dbsnp+' -O ' + resdir + '/' + pollo + '_recal_data.table1'
	cdata['dependency'] = 'afterok:' + str(p)
	p = send_sbatch(cdata)
	# ApplyBQSR, depende de BaseRecalibrator
	cdata['job_name'] = pollo + '_ApplyBQSR'
	cdata['filename'] = outdir + '/' + pollo + '_ApplyBQSR.sh'
	cdata['output'] = outdir + '/' + pollo + '_ApplyBQSR.out'
	cdata['command'] = gatk + ' ApplyBQSR -R ' + ref_fa+' -I ' + resdir + '/' + pollo + '_rmdups.bam -bqsr-recal-file ' + resdir + '/' + pollo + '_recal_data.table1 -O ' + resdir + '/' + pollo + '_recal.bam'
	cdata['dependency'] = 'afterok:' + str(p)
	p0 = send_sbatch(cdata)
	# AnalyzeCovariates, depende de BaseRecalibrator
	cdata['job_name'] = pollo + '_AnalyzeCovariates_1'
	cdata['filename'] = outdir + '/' + pollo + '_AnalyzeCovariates_1.sh'
	cdata['output'] = outdir + '/' + pollo + '_AnalyzeCovariates_1.out'
	cdata['command'] = gatk + ' AnalyzeCovariates -bqsr ' + resdir + '/' + pollo + '_recal_data.table1 --plots ' + resdir + '/' + pollo + '_AnalyzeCovariates.pdf'
	cdata['dependency'] = 'afterok:' + str(p)
	p1 = send_sbatch(cdata)
	ljobids.append(p1)
	# BaseRecalibrator 2, depende de BaseRecalibrator
	cdata['job_name'] = pollo + '_BaseRecalibrator_2'
	cdata['filename'] = outdir + '/' + pollo + '_BaseRecalibrator_2.sh'
	cdata['output'] = outdir + '/' + pollo + '_BaseRecalibrator_2.out'
	cdata['command'] = gatk + ' BaseRecalibrator -I '+ resdir + '/' + pollo + '_recal.bam -R '+ ref_fa+' --known-sites ' + ref_dir + '/'+known1+' --known-sites '+ ref_dir + '/'+known2+' --known-sites ' + ref_dir + '/'+dbsnp+' -O ' + resdir + '/' + pollo + '_recal_data.table2'
	cdata['dependency'] = 'afterok:' + str(p0)
	p1 = send_sbatch(cdata)
	# AnalyzeCovariates 2, depende de BaseRecalibrator 2
	cdata['job_name'] = pollo + '_AnalyzeCovariates_2'
	cdata['filename'] = outdir + '/' + pollo + '_AnalyzeCovariates_2.sh'
	cdata['output'] = outdir + '/' + pollo + '_AnalyzeCovariates_2.out'
	cdata['command'] = gatk + ' AnalyzeCovariates -before '+  resdir + '/' + pollo + '_recal_data.table1 -after '+  resdir + '/' + pollo + '_recal_data.table2 -plots '+  resdir + '/' + pollo +'_before-after-plots.pdf'
	cdata['dependency'] = 'afterok:' + str(p1)
	p1 = send_sbatch(cdata)
	ljobids.append(p1)
	# CollectHSMetrics 2, depende de ApplyBQSR
	cdata['job_name'] = pollo + '_CollectHsMetrics_2'
	cdata['filename'] = outdir + '/' + pollo + '_CollectHsMetrics_2.sh'
	cdata['output'] = outdir + '/' + pollo + '_CollectHsMetrics_2.out'
	cdata['command'] = gatk + ' CollectHsMetrics -BI '+baits + ' -TI ' + targets + ' -I ' + resdir + '/' + pollo + '_recal.bam -O ' + resdir + '/' + pollo + '_recal_hsmetrics.txt'
	cdata['dependency'] = 'afterok:' + str(p0)
	p1 = send_sbatch(cdata)
	ljobids.append(p1)
	# HaplotypeCaller, depende de ApplyBQSR
	cdata['job_name'] = pollo + '_HaplotypeCaller'
	cdata['filename'] = outdir + '/' + pollo + '_HaplotypeCaller.sh'
	cdata['output'] = outdir + '/' + pollo + '_HaplotypeCaller.out'
	cdata['command'] = gatk + ' HaplotypeCaller -R '+ref_fa+' -I '+ resdir + '/' + pollo + '_recal.bam -ERC GVCF --dbsnp '+ref_dir+'/'+dbsnp+' -O '+ resdir + '/' + pollo +'_raw.snps.indels.g.vcf.gz\n'
	cdata['command'] += gatk + ' VariantEval -R '+ref_fa+' -L '+unions+' -D '+ref_dir+'/'+hcsnps+' -O '+resdir + '/' + pollo +'_eval.gatkreport --eval '+ resdir + '/' + pollo +'_raw.snps.indels.g.vcf.gz'
	cdata['dependency'] = 'afterok:' + str(p0)
	p1 = send_sbatch(cdata)
	ljobids.append(p1)
	cdata['job_name'] = pollo + '_close_sbj'
	cdata['filename'] = outdir + '/' + pollo + '_close_sbj.sh'
	cdata['output'] = outdir + '/' + pollo + '_close_sbj.out'
	sdeps = 'afterok:'+',afterok:'.join(map(str,ljobids))
	cdata['dependency'] = sdeps
	cdata['command'] = 'rm -rf '+tmpdir
	sp = send_sbatch(cdata)
	jobids.append(sp)
deps = 'afterok:'+',afterok:'.join(map(str,jobids))
wdata = {'job_name':'end_wes', 'filename':outdir+'/warning_end.sh','output':outdir+'/warning_end.out', 'dependency':deps}
send_sbatch(wdata)
