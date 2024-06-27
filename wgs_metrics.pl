#!/usr/bin/perl 
use strict;
use warnings;
use SLURMACE qw(send2slurm);
use Data::Dump qw(dump);

################### Paths ##################
my $ipath = '/nas/Genomica/01-Data/02-WXS/02-Processed/202211_WES_PSP-DEGESCO/output/';
my $intervals = '/nas/Genomica/01-Data/00-Reference_files/02-GRCh38/202211_KAPA_WES_PSP-DEGESCO/KAPA_HyperExome_hg38_union.interval_list';
my $href = '/nas/Genomica/01-Data/00-Reference_files/02-GRCh38/00_Bundle/Homo_sapiens_assembly38.fasta';
my $gatk = 'singularity run --cleanenv -B /nas:/nas -B /ruby:/ruby -B /the_dysk:/the_dysk /nas/usr/local/opt/gatk4.simg gatk --java-options "-Djava.io.tmpdir='.$ENV{'TMPDIR'}.' -DGATK_STACKTRACE_ON_USER_EXCEPTION=true -Xmx16G"';
my $wdir = 'slurm';

############ Busca las muestras  #########################
my @idirs = glob( $ipath.'*' );
my @pollos;
foreach my $idir (@idirs) {
	if ( -d $idir ) {
		my ($pollo) = $idir =~ /.*\/output\/(.*)$/;
		push @pollos, $pollo;
	}
}

my %ptask = ( 'job_name' => 'wgs_metrics',
	'cpus' => 4,
	'time' => '24:0:0',
	'mailtype' => 'FAIL,TIME_LIMIT,STAGE_OUT',
	'partition' =>'fast',
);

foreach my $pollo (@pollos) {
	$ptask{'filename'} = $wdir.'/'.$pollo.'.sh';
	$ptask{'output'} = $wdir.'/'.$pollo.'-%j';
	$ptask{'command'} = $gatk.' CollectWgsMetrics --INPUT '.$ipath.$pollo.'/results/'.$pollo.'_recal.bam  --OUTPUT '.$ipath.$pollo.'/results/'.$pollo.'_wgs_metrics.txt --REFERENCE_SEQUENCE '.$href.' --COUNT_UNPAIRED true -Q 20 -MQ 0 --INTERVALS '.$intervals."\n";
	$ptask{'command'} .= $gatk.' CollectWgsMetrics --INPUT '.$ipath.$pollo.'/results/'.$pollo.'_recal.bam  --OUTPUT '.$ipath.$pollo.'/results/'.$pollo.'_raw_wgs_metrics.txt --REFERENCE_SEQUENCE '.$href.' --COUNT_UNPAIRED true -Q 3 -MQ 0 --INTERVALS '.$intervals;
	send2slurm(\%ptask);
}

my %warn = ('job_name' => 'wgs_metrics',         
	'filename' => $wdir.'/tasks_end.sh',         
	'mailtype' => 'END',         
	'output' => $wdir.'/tasks_end',         
	'dependency' => 'singleton' 
); 
send2slurm(\%warn);
