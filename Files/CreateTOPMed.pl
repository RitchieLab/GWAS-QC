#!/usr/bin/perl
#
# Program to reformat the TOPMed vcf to a format for pre-imputation checking.
# W.Rayner 2020
# william.rayner@helmholtz-muenchen.de
# wrayner@well.ox.ac.uk
# 
# v1.0.0 Created
#
#
#
#
#
#
#
#
#

use strict;
use warnings;
use File::Basename;
use File::Spec;
use Getopt::Long;
use Cwd;
use IO::Uncompress::Gunzip qw(gunzip $GunzipError);
use IO::Compress::Gzip qw(gzip $GzipError);

my $in;
my $out;
my $z;
my $f = 0;
my $all = 0;

GetOptions
 (
 "i|infile=s"      => \$in,      # vcf file containing the TOPMed frequencies
 "o|outfile=s"     => \$out,     # Set the output filename to override the default
 "a|all"           => \$all      # flag to output all variants not just PASS variants
 );

my ($volume,$directories,$parsed_file) = File::Spec->splitpath($in);
my $abs_path = File::Spec->rel2abs($directories);

if (!$in)
 {
 print "ERROR: No input file specified\n";
 exit;
 }

if (!-e $in)
 {
 print "ERROR: Can't find file specified:\n $in\n";
 exit;
 }

if (!$out)
 {
 #create output file
 if ($parsed_file =~ /(.*?)\.vcf.*/)
  {
  my $fileStem = $1;
  if (!$all)
   {
   $fileStem =~ s/ALL//g;
   $fileStem = 'PASS.Variants'.$fileStem;
   }
  $out = $abs_path.'/'.$fileStem.'.tab.gz';
  }
 }

print "Input File:  $in\nOutput File: $out\n";
#open O, ">$out" or die $!;
my $zo = new IO::Compress::Gzip "$out" or die "IO::Compress::Gzip failed: $GzipError\n"; 
print $zo "#CHROM\tPOS\tID\tREF\tALT\tAC\tAN\tAF\n";

$parsed_file = $abs_path.'/'.$parsed_file;
print "$parsed_file\n";

if ($parsed_file =~ /.*\.gz$/)
 {
 if (-e '/bin/gunzip')
  {
  open $z, '-|', '/bin/gunzip', '-c', "$parsed_file";
  }
 else
  {
  print "Gunzip not located, trying alternative method which may fail on bgzipped files\n";
  $z = new IO::Uncompress::Gunzip "$parsed_file" or die "IO::Uncompress::Gunzip failed: $GunzipError\n"; 
  }
 }
else
 {
 open $z, "$parsed_file" or die $!;
 }
 
while (<$z>)
 {
 chomp;

 if ($f) 
  {
  my %data;
  s/^chr//;
  my @temp = split/\t/;
  my $filter = $temp[6];

  if ($all) #if exporting all then set filter to PASS for all variants
   {
   $filter = 'PASS';
   }

  if ($filter eq 'PASS')
   {
   my @entries = split(/\;/, $temp[7]);
   foreach my $entry (@entries)
    {
    $entry =~ /(.*)\=(.*)/;
    $data{$1} = $2;
    }
   
   print $zo "$temp[0]\t$temp[1]\t$temp[2]\t$temp[3]\t$temp[4]";
   if ($data{'AC'})
    {
    print $zo "\t$data{'AC'}";
    }
   else
    {
    print $zo "\t0";
    }

   if ($data{'AN'})
    {
    print $zo "\t$data{'AN'}";
    }
   else
    {
    print $zo "\t0";
    }

   if ($data{'AF'}) 
    {
    print $zo "\t$data{'AF'}";
    } 
   else
    {
    print $zo "\t0";
    }
   print $zo "\n";
   }
  }

 if (/#CHROM/)
  {
  $f = 1;
  }

 }

