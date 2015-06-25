#!/usr/bin/env perl
#
# Example code for generating a minimal VCF file using the perl API
#
# Author: pd3@sanger
#

use strict;
use warnings;
use Carp;
use Vcf;

my $sample = 'Sample1';
my $vcf_out = Vcf->new();
$vcf_out->add_columns($sample);
$vcf_out->add_header_line({key=>'FORMAT',ID=>'GT',Number=>'1',Type=>'String',Description=>"Genotype"});
$vcf_out->add_header_line({key=>'ALT',ID=>'DEL',Description=>"Deletion"});
$vcf_out->add_header_line({key=>'ALT',ID=>'DEL:ME:ALU',Description=>"Deletion of ALU element"});
$vcf_out->add_header_line({key=>'ALT',ID=>'DEL:ME:L1',Description=>"Deletion of L1 element"});
$vcf_out->add_header_line({key=>'ALT',ID=>'DUP',Description=>"Duplication"});
$vcf_out->add_header_line({key=>'INFO',ID=>'DP',Number=>1,Type=>'Integer',Description=>"Total Depth"});
$vcf_out->add_header_line({key=>'INFO',ID=>'H2',Number=>0,Type=>'Flag',Description=>"HapMap2 membership"});

print $vcf_out->format_header();

my $pos = 1;
for my $gt qw(A/A C/C <DEL>/C <DEL:ME:ALU>/<DEL:ME:ALU> <DEL:ME:L1>/<DEL:ME:L1> <DUP>/<DUP>)
{
    $pos++;

    my %out;
    $out{CHROM}  = '1';
    $out{POS}    = $pos;
    $out{ID}     = '.';
    $out{ALT}    = [];
    $out{REF}    = 'C';
    $out{QUAL}   = '.';
    $out{FILTER} = ['.'];
    $out{INFO}   = { DP=>3, H2=>undef };
    $out{FORMAT} = ['GT'];
    $out{gtypes}{$sample}{GT} = $gt;

    $vcf_out->format_genotype_strings(\%out);
    print $vcf_out->format_line(\%out);
}

