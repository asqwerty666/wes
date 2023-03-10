#!/usr/bin/perl
use strict;
use warnings;
use Data::Dump qw(dump);
use File::Find::Rule;
use File::Basename qw(basename);
my $tmp_dir = $ENV{'TMPDIR'};
my $wdir = '/ruby/WES_full_output';
my %rep;
my @vtags = ('Per base sequence quality',
	'Per sequence quality scores',
	'Per base sequence content',
	'Per sequence GC content',
	'Per base N content',
	'Sequence Length Distribution',
	'Sequence Duplication Levels',
	'Overrepresented sequences',
	'Adapter Content'
);
my @stags = ('Total Sequences',
	'Sequences flagged as poor quality',
	'Sequence length',
	'%GC'
);
my @sbj_file1 = find(file => 'name' => "*_1_fastqc.html", in => $wdir);
foreach my $sbjf (@sbj_file1){
	my $sbj = basename $sbjf;
	$sbj =~ s/(.*)_1_fastqc.html/$1/;
	foreach my $i (1,2){
		my $ifile = $wdir.'/'.$sbj.'/qc/'.$sbj.'_'.$i.'_fastqc.html';
		open IDF, "lynx -dump $ifile |";
		while(<IDF>){
			 foreach my $tag (@vtags){
				 if(/^\[[A-Z]+\]\s+$tag\s*$/){
					 my ($vcap) = /^\[([A-Z]+)\]\s+$tag\s*$/;
					 $rep{$sbj}{$i}{$tag} = $vcap;
				 }
			 }
			 foreach my $tag (@stags){
				 if(/^\s+$tag\s+\d+\s*$/){
					  my ($vcap) = /^\s+$tag\s+(\d+)\s*$/;
					  $rep{$sbj}{$i}{$tag} = $vcap;
				 }
			 }
		 }
		 close IDF;
	 }
}
print "FASTQC";
foreach my $tag (sort @stags){
	print ",$tag";
}
foreach my $tag (sort @vtags){
	print ",$tag"; 
}
print "\n";
foreach my $sbj (sort keys %rep){
	foreach my $chunk (sort keys %{$rep{$sbj}}){
		my $acc = $sbj.'_'.$chunk;
		print "$acc";
		foreach my $tag (sort @stags){
			if(exists($rep{$sbj}{$chunk}{$tag})){
				print ",$rep{$sbj}{$chunk}{$tag}";
			}else{
				print ",NULL";
			}
		}
		foreach my $tag (sort @vtags){
			if(exists($rep{$sbj}{$chunk}{$tag}) and $rep{$sbj}{$chunk}{$tag}){
				print ",$rep{$sbj}{$chunk}{$tag}";
			}else{
				print ",NULL";
			}
		}
		print "\n";
	}
}


