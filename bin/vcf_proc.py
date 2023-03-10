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
pollos_list = wdir+'/subjects.list'
if os.path.exists(pollos_list):  os.remove(pollos_list)
for pollo in dir_cont:
	fq_list = []
	fq_id =''
	fq_files = next(os.walk(lpath+'/'+pollo))[2]
	udir = wdir+'/'+pollo
	tmpdir = udir+'/tmp'
	resdir = udir+'/results'
	with open(pollos_list, 'a') as lf:
		lf.write(pollo+'\t'+resdir + '/' + pollo +'_raw.snps.indels.g.vcf.gz\n')
ldata = {'time':'72:0:0', 'cpus':'16', 'mem-per-cpu':'4G'}
chjobids = []
for ch in range(22):
	ldata['job_name'] = 'GenoypeGVCFs_'+str(ch+1)
	ldata['filename'] = outdir + '/' + 'GenoypeGVCFs_'+str(ch+1)+'.sh'
	ldata['output'] = outdir + '/' + 'GenoypeGVCFs_'+str(ch+1)+'.out'
#	ldata['dependency'] = deps
	dbdir = fwdir+'/wes.chr'+str(ch+1)+'.db'
	ldata['command'] = gatk +' GenomicsDBImport --genomicsdb-workspace-path '+dbdir+' --batch-size 50 --sample-name-map '+pollos_list+' -L chr'+str(ch+1)+' --reader-threads 8\n'
	ldata['command'] += gatk +' GenotypeGVCFs -R '+ref_fa+' -V gendb://'+dbdir+' -O '+fwdir+'/wes_chr'+str(ch+1)+'.snps.indels.g.vcf.gz\n'
	ldata['command'] += 'gunzip -c '+ fwdir+'/wes_chr'+str(ch+1)+'.snps.indels.g.vcf.gz |'+snpEff+' ann -chr '+str(ch+1)+' -s snpEff_summary_chr'+str(ch+1)+'.html hg38 - | gzip -c > '+fwdir+'/wes_chr'+str(ch+1)+'.snps.indels.g.ann.vcf.gz'
	chp = send_sbatch(ldata)
	chjobids.append(chp)
deps = 'afterok:'+',afterok:'.join(map(str,chjobids))
wdata = {'job_name':'end_wes', 'filename':outdir+'/warning_end.sh','output':outdir+'/warning_end.out', 'dependency':deps}
send_sbatch(wdata)
