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
src_dir = '/nas/GRACE/seqEOAD_DEGESCO/input/DEGESCO/muestras'
ref_dir = '/ruby/WES_reference/hg38'
ref_name = 'Homo_sapiens_assembly38'
ref_fa = ref_dir+'/'+ref_name+'.fasta'
output_dir = '/nas/osotolongo/wes1/haplos'
tmp_dir = '/ruby/'+os.environ.get('USER')+'/tmp/'
#baits = ref_dir + '/KAPA_HyperExome_hg38_bait.interval_list'
#targets = ref_dir + '/KAPA_HyperExome_hg38_target.interval_list'
#unions = ref_dir + '/KAPA_HyperExome_hg38_union.interval_list' 
#known1 = 'Homo_sapiens_assembly38.known_indels.vcf.gz'
#known2 = 'Mills_and_1000G_gold_standard.indels.hg38.vcf.gz'
dbsnp = 'dbsnp_138.b37.vcf'
#hcsnps = '1000G_phase1.snps.high_confidence.hg38.vcf.gz'
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
hwdir =  wdir+'/reverts'
if not os.path.isdir(hwdir): os.mkdir(hwdir)
fq = {}
jobids = []
cdata = {'time':'24:0:0', 'cpus':'4', 'mem_cpu':'4G'}
cdata['test'] = 1
pollos_list = wdir+'/subjects.list'
if os.path.exists(pollos_list):  os.remove(pollos_list)
for pollo in dir_cont:
    tmp_dir = tmp_shit+'/'+pollo
    #if not os.path.isdir(tmp_dir): os.mkdir(tmp_dir)
    ljobids = []
    call = [samtools, 'view', '-H', src_dir + '/' + pollo + '/' + pollo +'.bam']
    result = subprocess.run(call, stdout=subprocess.PIPE)
    rdata = list(filter(lambda x:'@RG' in x, str(result.stdout).split('\\n')))[0].split('\\t')
    my_id=list(filter(lambda x:'ID:' in x, rdata))[0].split(':')[1]
    with open(pollos_list, 'a') as lf:
        lf.write(pollo+'\t'+hwdir + '/' + pollo + '.snps.indels.g.vcf.gz\n')
    cdata['job_name'] = pollo + '_RevertSam'
    cdata['filename'] = outdir + '/' + pollo + '_RevertSam.sh'
    cdata['output'] = outdir + '/' + pollo + '_RevertSam.out'
    cdata['command'] = 'mkdir '+tmp_dir+' \n'
    cdata['command'] += 'cp ' + src_dir + '/' + pollo + '/' + pollo +'.b* '+ tmp_dir +'/\n'
    cdata['command'] += gatk + ' RevertSam -I '+ tmp_dir + '/' + pollo +'.bam -O '+tmp_dir + '/' + pollo +'_u.bam\n'
    cdata['command'] += gatk + ' SortSam -I ' + tmp_dir + '/' + pollo +'_u.bam -O ' + tmp_dir + '/' + pollo +'_us.bam -SORT_ORDER queryname\n'
    cdata['command'] += gatk + ' SamToFastq -I ' + tmp_dir + '/' + pollo +'.bam -FASTQ ' + tmp_dir + '/' + pollo + '_a.fastq -F2 ' +  tmp_dir + '/' + pollo + '_b.fastq\n'
    cdata['command'] += bwa + ' -R "@RG\\tID:'+my_id+'\\tPL:ILLUMINA\\tLB:CUSTOM\\tSM:'+pollo+'" ' + ref_fa + ' ' + tmp_dir + '/' +  pollo + '_a.fastq ' + tmp_dir + '/' +  pollo + '_b.fastq > ' + tmp_dir + '/' + pollo +'_alignment.sam\n'
    cdata['command'] += gatk + ' SortSam -I ' + tmp_dir + '/' + pollo +'_alignment.sam -O ' + tmp_dir + '/' + pollo +'_salignment.sam -SORT_ORDER queryname\n'
    cdata['command'] += gatk + ' MergeBamAlignment -ALIGNED ' + tmp_dir + '/' + pollo +'_salignment.sam -UNMAPPED ' + tmp_dir + '/' + pollo +'_us.bam -R ' + ref_fa + ' -O ' + hwdir + '/' + pollo + '.sam\n'
    cdata['command'] += 'rm -rf '+tmp_dir
    #cdata['dependency'] = 'afterok:' + str(p0)
    p1 = send_sbatch(cdata)
    ljobids.append(p1)
deps = 'afterok:'+',afterok:'.join(map(str,ljobids))
if not 'test' in cdata or not cdata['test']:
    wdata = {'job_name':'end_wes', 'filename':outdir+'/warning_end.sh','output':outdir+'/warning_end.out', 'dependency':deps}
    send_sbatch(wdata)
