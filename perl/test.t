#!/usr/bin/env perl
#
# Author: petr.danecek@sanger
#
# Usage: test.t [-d]
#

use strict;
use warnings;
use Carp;
use IPC::Open2;
use FindBin;
use lib "$FindBin::Bin";
use Vcf;

BEGIN {
    use Test::Most tests => 75;
}


my $path = $FindBin::RealBin;

my $debug = ($ARGV[0] && $ARGV[0] eq '-d') ? 1 : 0;

test_bgzip_and_tabix("$path/../examples/merge-test-a.vcf");
test_validator($path,"$path/../examples/valid-3.3.vcf");
test_validator($path,"$path/../examples/valid-4.0.vcf");
test_validator($path,"$path/../examples/valid-4.1.vcf");
test_validator($path,"$path/../examples/floats.vcf");
test_format_validation($path,'3.3');
test_format_validation($path,'4.0');
test_format_validation($path,'4.1');
test_parse($path);
test_vcf_stats($path,"$path/../examples/valid-4.0.vcf");
test_empty_cols($path,'4.0');
test_merge($path,'merge-test.vcf.out','merge-test-a.vcf','merge-test-b.vcf','merge-test-c.vcf');
test_compare($path,'cmp-test-a.vcf','cmp-test-b.vcf','cmp-test.out');
test_isec($path,'-n +2','isec-n2-test.vcf.out','merge-test-a.vcf','merge-test-b.vcf','merge-test-c.vcf');
test_query_vcf("$path/../examples/",'cmp-test-a.vcf','query-test.out','%CHROM:%POS\tref=%REF\talt=%ALT\tqual=%QUAL\t%INFO/DP[\t%SAMPLE=%GT]\n');
test_shuffle("$path/../examples/",'cmp-test-a.vcf','shuffle-test.vcf');
test_concat("$path/../examples/",'concat.out','concat-a.vcf','concat-b.vcf','concat-c.vcf');
test_annotate("$path/../examples/",'-c FROM,TO,CHROM,-,-,-,INFO/HM2,INFO/GN,INFO/DP -d key=INFO,ID=HM2,Number=0,Type=Flag,Description="HapMap2 membership" -d key=INFO,ID=GN,Number=1,Type=String,Description="Gene Name" -d key=INFO,ID=DP,Number=0,Type=Integer,Description="Depth,etc"','annotate.out','concat-a.vcf','annotate.txt');
test_annotate("$path/../examples/",'-c FROM,TO,CHROM,ID,REF,ALT,INFO/HM2,INFO/GN,INFO/DP -d key=INFO,ID=HM2,Number=0,Type=Flag,Description="HapMap2 membership" -d key=INFO,ID=GN,Number=1,Type=String,Description="Gene Name" -d key=INFO,ID=DP,Number=0,Type=Integer,Description="Depth,etc"','annotate3.out','concat-a.vcf','annotate.txt');
test_annotate("$path/../examples/",'-f +/D=34/c=2,3','annotate2.out','annotate-test.vcf');
test_fill_an_ac("$path/../examples/",'fill-an-ac.out','concat-a.vcf');
test_indel_stats("$path/../examples/",'indel-stats.out','indel-stats.vcf','indel-stats.tab');
test_consensus("$path/../examples/",'','consensus.out','consensus.vcf','consensus.fa');
test_consensus("$path/../examples/",'-s NA001','consensus.out2','consensus.vcf','consensus.fa');
test_contrast("$path/../examples/",'-n +D -A,B,C -d 10','contrast.out','contrast.vcf');
test_ploidy("$path/../examples/",'fix-ploidy');
test_api_event_type([qw(A C),'s 1 C'],[qw(A ACGT),'i 3 CGT'],[qw(ACGT A),'i -3 CGT'],[qw(ACGT ACT),'i -1 G'],
    [qw(ACGT AAA),'o 3 AAA'],[qw(A .),'r 0 A'],[qw(A <ID>),'u 0 <ID>'],[qw(ACG AGC),'s 2 AGC'], [qw(A .A),'b'], [qw(A A.),'b']);
test_api();

exit;

#--------------------------------------

sub test_bgzip_and_tabix
{
    my ($file) = @_;
    my $cmd;

    $cmd = "cat $file | bgzip -c > $file.gz";
    system($cmd);
    is($?,0,"Is bgzip OK? .. $cmd");

    $cmd = "tabix $file.gz";
    system($cmd);
    is($?,0,"Is tabix OK? .. $cmd");
}

sub test_validator
{
    my ($path,$fname) = @_;
    
    my $cmd = "perl -I$path -MVcf -e validate $fname";
    my @out = `$cmd 2>&1`;
    my @exp = ();
    is_deeply(\@out,\@exp,"Testing validator .. $cmd");
}

sub test_format_validation
{
    my ($path,$version) = @_;

    my ($chld_in,$chld_out);
    my $cmd = "perl -I$path -MVcf -e validate 2>&1";
    my $pid = open2($chld_out, $chld_in, $cmd);

    my $vcf = Vcf->new(version=>$version);
    $vcf->recalc_ac_an(2);
    $vcf->add_header_line({key=>'INFO', ID=>'AC',Number=>-1,Type=>'Integer',Description=>'Allele count in genotypes'});
    $vcf->add_header_line({key=>'INFO', ID=>'AN',Number=>1,Type=>'Integer',Description=>'Total number of alleles in called genotypes'});
    $vcf->add_header_line({key=>'FORMAT', ID=>'GT',Number=>1,Type=>'String',Description=>'Genotype'});
    if ( $version >= 4.0 )
    {
        $vcf->add_header_line({key=>'ALT',ID=>'DEL:ME:ALU', Description=>'Deletion of ALU element'});
    }
    if ( $version >= 4.1 )
    {
        $vcf->add_header_line({key=>'reference',value=>'file:/some/file.fa'});
        $vcf->add_header_line({key=>'contig',ID=>'1',length=>12345,md5=>'f126cdf8a6e0c7f379d618ff66beb2da',assembly=>'E.T.'});
    }
    $vcf->add_columns('NA0001','NA0002');
    print $vcf->format_header() unless !$debug;
    print $chld_in $vcf->format_header();

    my %rec = ( CHROM=>1, POS=>1, REF=>'A', QUAL=>$$vcf{defaults}{QUAL}, FORMAT=>['GT'] );
    $rec{gtypes}{NA0001}{GT} = 'A/A';
    $rec{gtypes}{NA0002}{GT} = $$vcf{defaults}{GT};
    $vcf->format_genotype_strings(\%rec);
    print $vcf->format_line(\%rec) unless !$debug;
    print $chld_in $vcf->format_line(\%rec);

    $rec{POS} = 2;
    $rec{gtypes}{NA0002}{GT} = 'IA|D1';
    if ( $version >= 4.0 ) 
    { 
        $rec{REF} = 'AC';
        $rec{gtypes}{NA0002}{GT} = 'ATC|<DEL:ME:ALU>'; 
    }
    $vcf->format_genotype_strings(\%rec);
    print $vcf->format_line(\%rec) unless !$debug;
    print $chld_in $vcf->format_line(\%rec);
    close($chld_in);

    my @exp = ();
    my @out = ();
    while (my $line=<$chld_out>)
    {
        chomp($line);
        push @out,$line;
    }
    close($chld_out);
    waitpid $pid, 0;

    if ( !is_deeply(\@out,\@exp,"Testing formatting followed by validation .. $cmd") )
    {
        print STDERR @out;
    }
}

sub test_parse
{
    my ($path) = @_;
    my $vcf = Vcf->new(file=>"$path/../examples/parse-test.vcf");
    $vcf->parse_header;
    my $line;
    $line = $vcf->next_data_array; is_deeply($$line[4],"G","Testing next_data_array");
    $line = $vcf->next_data_array; is_deeply($$line[4],"G,<DEL2>,T,<DEL3>","Testing next_data_array");
    $line = $vcf->next_data_array; is_deeply($$line[4],"<DEL1>,G,<DEL2>,T","Testing next_data_array");
    $line = $vcf->next_data_array; is_deeply($$line[4],"<DEL1>,G,<DEL2>,T,<DEL3>","Testing next_data_array");
}

sub test_vcf_stats
{
    my ($path,$file) = @_;
    my $cmd = "perl -I$path -MVcf $path/vcf-stats $file";
    my @out = `$cmd 2>&1`;
    open(my $fh,'<',"$file.stats") or confess("$file.stats: $!");
    my @exp = <$fh>;
    close($fh);

    is_deeply(\@out,\@exp,"Testing vcf-stats .. $cmd");
}

sub test_empty_cols
{
    my ($path,$version) = @_;

    my ($header,$vcf,@out,$exp);

    $vcf = Vcf->new(version=>$version);
    $vcf->add_header_line({key=>'FORMAT', ID=>'GT',Number=>1,Type=>'String',Description=>'Genotype'});
    $vcf->add_columns(qw(CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  NA0001));
    $header = $vcf->format_header();
    @out = split(/\n/,$header);
    $exp = join("\t",qw(CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  NA0001));
    is_deeply($out[-1],'#'.$exp,"Testing add_columns with genotypes full, $version.");

    $vcf = Vcf->new(version=>$version);
    $vcf->add_header_line({key=>'FORMAT', ID=>'GT',Number=>1,Type=>'String',Description=>'Genotype'});
    $vcf->add_columns('NA0001');
    $header = $vcf->format_header();
    @out = split(/\n/,$header);
    $exp = join("\t",qw(CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  NA0001));
    is_deeply($out[-1],'#'.$exp,"Testing add_columns with genotypes brief, $version.");

    $vcf = Vcf->new(version=>$version);
    $vcf->add_header_line({key=>'FORMAT', ID=>'GT',Number=>1,Type=>'String',Description=>'Genotype'});
    $vcf->add_columns();
    $header = $vcf->format_header();
    @out = split(/\n/,$header);
    $exp = join("\t",qw(CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO));
    is_deeply($out[-1],'#'.$exp,"Testing add_columns brief, $version.");

    $vcf = Vcf->new(version=>$version);
    $vcf->add_header_line({key=>'FORMAT', ID=>'GT',Number=>1,Type=>'String',Description=>'Genotype'});
    $vcf->add_columns('FORMAT');
    $header = $vcf->format_header();
    @out = split(/\n/,$header);
    $exp = join("\t",qw(CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO FORMAT));
    is_deeply($out[-1],'#'.$exp,"Testing add_columns no gtypes, $version.");
}

sub test_compare
{
    my ($path,$a,$b,$expected) = @_;

    my $curdir = `pwd`;
    chomp($curdir);
    chdir("$path/../examples");

    for my $file ($a,$b)
    {
        `cat $file | bgzip -c > $file.gz`;
        `tabix -p vcf -f $file.gz`;
    }
    
    my $cmd = "perl -I../perl/ -MVcf ../perl/vcf-compare -g $a.gz $b.gz | grep -v '^# The command'";
    my @out = `$cmd 2>&1`;
    open(my $fh,'<',"$expected") or confess("$expected: $!");
    my @exp = <$fh>;
    close($fh);

    chdir($curdir);

    is_deeply(\@out,\@exp,"Testing vcf-compare .. $cmd");
}

sub test_merge
{
    my ($path,$expected,@files) = @_;

    my $curdir = `pwd`;
    chomp($curdir);
    chdir("$path/../examples");

    my $cmd = "perl -I../perl/ -MVcf ../perl/vcf-merge";
    for my $file (@files)
    {
        `cat $file | bgzip -c > $file.gz; tabix -f -p vcf $file.gz`;
        $cmd .= " $file.gz";
    }
    my @out = `$cmd 2>/dev/null | grep -v ^##source`;
    open(my $fh,'<',$expected) or confess("$expected: $!");
    my @exp = <$fh>;
    close($fh);

    chdir($curdir);
    is_deeply(\@out,\@exp,"Testing vcf-merge .. $cmd");
}

sub test_isec
{
    my ($path,$opts,$expected,@files) = @_;

    my $curdir = `pwd`;
    chomp($curdir);
    chdir("$path/../examples");

    my $cmd = "perl -I../perl/ -MVcf ../perl/vcf-isec -f $opts";
    for my $file (@files)
    {
        `cat $file | bgzip -c > $file.gz; tabix -f -p vcf $file.gz`;
        $cmd .= " $file.gz";
    }
    my @out = `$cmd 2>&1 | grep -v ^##source`;
    open(my $fh,'<',$expected) or confess("$expected: $!");
    my @exp = <$fh>;
    close($fh);

    chdir($curdir);
    is_deeply(\@out,\@exp,"Testing vcf-isec .. $cmd");
}


sub test_query_vcf
{
    my ($path,$file,$expected,$query) = @_;

    my $curdir = `pwd`;
    chomp($curdir);
    chdir("$path/../examples");

    my $cmd = "perl -I../perl/ -MVcf ../perl/vcf-query -f '$query' $file";
    my @out = `$cmd 2>&1`;
    open(my $fh,'<',$expected) or confess("$expected: $!");
    my @exp = <$fh>;
    close($fh);

    chdir($curdir);
    is_deeply(\@out,\@exp,"Testing vcf-query .. $cmd");
}


sub test_shuffle
{
    my ($path,$template,$file) = @_;

    my $curdir = `pwd`;
    chomp($curdir);
    chdir("$path/../examples");

    my $cmd = "perl -I../perl/ -MVcf ../perl/vcf-shuffle-cols -t $template $file";
    my @out = `$cmd 2>&1`;
    open(my $fh,'<',$template) or confess("$template: $!");
    my @exp = <$fh>;
    close($fh);

    chdir($curdir);
    is_deeply(\@out,\@exp,"Testing vcf-shuffle-cols .. $cmd");
}

sub test_concat
{
    my ($path,$out,@files) = @_;

    my $curdir = `pwd`;
    chomp($curdir);
    chdir("$path/../examples");

    my $cmd = "perl -I../perl/ -MVcf ../perl/vcf-concat -s 3";
    for my $file (@files)
    {
        `cat $file | bgzip -c > $file.gz`;
        `tabix -p vcf -f $file.gz`;
        $cmd .= " $file.gz";
    }
 
    my @out = `$cmd 2>&1`;
    open(my $fh,'<',$out) or confess("$out: $!");
    my @exp = <$fh>;
    close($fh);

    chdir($curdir);
    is_deeply(\@out,\@exp,"Testing vcf-concat .. $cmd");
}


sub test_annotate
{
    my ($path,$args,$out,$vcf,$annot) = @_;

    my $curdir = `pwd`;
    chomp($curdir);
    chdir("$path/../examples");

    my $cmd = "perl -I../perl/ -MVcf ../perl/vcf-annotate $args $vcf";

    if ( defined $annot )
    {
        `cat $annot | bgzip -c > $annot.gz`;
        `tabix -s 3 -b 1 -e 2 -f $annot.gz`;
        $cmd .= " -a $annot.gz";
    }

    my @out = `$cmd 2>&1 | grep -v ^##source`;
    open(my $fh,'<',$out) or confess("$out: $!");
    my @exp = <$fh>;
    close($fh);

    chdir($curdir);
    is_deeply(\@out,\@exp,"Testing vcf-annotate .. $cmd");
}

sub test_fill_an_ac
{
    my ($path,$out,$vcf) = @_;

    my $curdir = `pwd`;
    chomp($curdir);
    chdir("$path/../examples");

    my $cmd = "perl -I../perl/ -MVcf ../perl/fill-an-ac $vcf";
    my @out = `$cmd 2>&1`;
    open(my $fh,'<',$out) or confess("$out: $!");
    my @exp = <$fh>;
    close($fh);

    chdir($curdir);
    is_deeply(\@out,\@exp,"Testing fill-an-ac .. $cmd");
}

sub test_indel_stats
{
    my ($path,$out,$vcf,$tab) = @_;

    my $curdir = `pwd`;
    chomp($curdir);
    chdir("$path/../examples");

    my $cmd = "perl -I../perl/ -MVcf ../perl/vcf-indel-stats -e $tab < $vcf";
    my @out = `$cmd 2>&1`;
    open(my $fh,'<',$out) or confess("$out: $!");
    my @exp = <$fh>;
    close($fh);

    chdir($curdir);
    is_deeply(\@out,\@exp,"Testing fill-an-ac .. $cmd");
}

sub test_consensus
{
    my ($path,$args,$out,$vcf,$fa) = @_;

    my $curdir = `pwd`;
    chomp($curdir);
    chdir("$path/../examples");
    `cat $vcf | bgzip -c > $vcf.gz`;
    `tabix -p vcf -f $vcf.gz`;
    my $cmd = "perl -I../perl/ -MVcf ../perl/vcf-consensus $args $vcf.gz < $fa";
    my @out = `$cmd`;
    open(my $fh,'<',$out) or confess("$out: $!");
    my @exp = <$fh>;
    close($fh);

    chdir($curdir);
    is_deeply(\@out,\@exp,"Testing vcf-consensus .. $cmd");
}

sub test_contrast
{
	my ($path,$args,$out,$vcf) = @_;
	my $curdir = `pwd`;
	chomp($curdir);
	chdir("$path/../examples");
	my $cmd = "perl -I../perl/ -MVcf ../perl/vcf-contrast $args $vcf | grep -v ^##source";
	my @out = `$cmd 2>&1`;
	open(my $fh,'<',$out) or confess("$out: $!");
	my @exp = <$fh>;
	close($fh);

	chdir($curdir);
	is_deeply(\@out,\@exp,"Testing vcf-contrast .. $cmd");
}

sub test_ploidy
{
    my ($path,$prefix) = @_;
    my $curdir = `pwd`;
    chomp($curdir);
    chdir("$path/../examples");
    my $cmd = "cat $prefix.vcf | perl -I../perl/ -MVcf ../perl/vcf-fix-ploidy -s $prefix.samples -p $prefix.txt 2>/dev/null | vcf-query -f '\%POS[\\t\%SAMPLE \%GTR \%PL]\\n'";
    my @out = `$cmd 2>&1`;
    open(my $fh,'<',"$prefix.out") or confess("$prefix.out: $!");
    my @exp = <$fh>;
    close($fh);

    chdir($curdir);
    is_deeply(\@out,\@exp,"Testing vcf-fix-ploidy .. $cmd");
}

sub test_api_event_type
{
    my (@subs) = @_;
    my $vcf = Vcf->new();
    for my $mut (@subs)
    {
        my $exp = join(' ', $vcf->event_type($$mut[0],$$mut[1]));
        is_deeply($$mut[2],$exp,"Testing API event_type($$mut[0],$$mut[1]) .. $exp");
    }
}

sub test_api
{
    my $vcf = Vcf->new();

    my $ret;
    my $fmt = 'GT:GL:PL'; 
    $ret = $vcf->get_tag_index($fmt,'GT',':'); is($ret,0,"Testing get_tag_index($fmt,'GT',':')");
    $ret = $vcf->get_tag_index($fmt,'GL',':'); is($ret,1,"Testing get_tag_index($fmt,'GL',':')");
    $ret = $vcf->get_tag_index($fmt,'PL',':'); is($ret,2,"Testing get_tag_index($fmt,'PL',':')");

    $ret = $vcf->remove_field($fmt,0,':'); is($ret,'GL:PL',"Testing get_tag_index($fmt,0,':')");
    $ret = $vcf->remove_field($fmt,1,':'); is($ret,'GT:PL',"Testing get_tag_index($fmt,1,':')");
    $ret = $vcf->remove_field($fmt,2,':'); is($ret,'GT:GL',"Testing get_tag_index($fmt,2,':')");

    $ret = $vcf->replace_field($fmt,'XX',0,':'); is($ret,'XX:GL:PL',"Testing get_tag_index($fmt,'XX',0,':')");
    $ret = $vcf->replace_field($fmt,'XX',1,':'); is($ret,'GT:XX:PL',"Testing get_tag_index($fmt,'XX',1,':')");
    $ret = $vcf->replace_field($fmt,'XX',2,':'); is($ret,'GT:GL:XX',"Testing get_tag_index($fmt,'XX',2,':')");
    $ret = $vcf->replace_field($fmt,'XX',4,':'); is($ret,'GT:GL:PL::XX',"Testing get_tag_index($fmt,'XX',4,':')");

    $ret = $vcf->decode_genotype('C',[qw(G T)],'0/1/2|1/0|1|2'); is($ret,'C/G/T|G/C|G|T',"Testing decode_genotype('C',['G','T'],'0/1/2|1/0|1|2')");
    $ret = $vcf->decode_genotype('C',[qw(G T)],'2|1'); is($ret,'T|G',"Testing decode_genotype('C',['G','T'],'2|1')");
    $ret = $vcf->decode_genotype('C',[qw(G T)],'2'); is($ret,'T',"Testing decode_genotype('C',['G','T'],'2')");

    my $info = 'NS=2;HM;AF=0.333;AFA=T;DB';
    $ret = $vcf->get_info_field($info,'NS');  is($ret,'2',"Testing get_info_field($info,'NS')");
    $ret = $vcf->get_info_field($info,'AF');  is($ret,'0.333',"Testing get_info_field($info,'AF')");
    $ret = $vcf->get_info_field($info,'AFA'); is($ret,'T',"Testing get_info_field($info,'AFA')");
    $ret = $vcf->get_info_field($info,'HM');  is($ret,'1',"Testing get_info_field($info,'HM')");
    $ret = $vcf->get_info_field($info,'DB');  is($ret,'1',"Testing get_info_field($info,'DB')");
    $ret = $vcf->get_info_field($info,'DBX'); is($ret,undef,"Testing get_info_field($info,'DBX')");
    $ret = $vcf->get_info_field('DB','DB'); is($ret,'1',"Testing get_info_field('DB','DB')");
    $ret = $vcf->get_info_field('XDB','DB'); is($ret,undef,"Testing get_info_field('XDB','DB')");

    my @ret;
    @ret = $vcf->split_gt('0/1'); is_deeply(\@ret,[0,1],"Testing split_gt('0/1')");
    @ret = $vcf->split_gt('0'); is_deeply(\@ret,[0],"Testing split_gt('0')");

    my @als;
    @als = ("TTGGTAT","TTGGTATCTAGTGGTAT,TGGTATCTAGTGGTAT"); @ret = $vcf->normalize_alleles(@als); 
    is_deeply(\@ret,["T","TTGGTATCTAG","TGGTATCTAG"],"Testing normalize_alleles(".join(',',@als).")");
    @als = ("TT","TCTAGTGGTAAT,TCT"); @ret = $vcf->normalize_alleles(@als); 
    is_deeply(\@ret,["T","TCTAGTGGTAA","TC"],"Testing normalize_alleles(".join(',',@als).")");
    @als = ("TGGGGGG","TGGGGGGG"); @ret = $vcf->normalize_alleles(@als); 
    is_deeply(\@ret,["T","TG"],"Testing normalize_alleles(".join(',',@als).")");
    @als = ("CAAAAAA","CAAAAA"); @ret = $vcf->normalize_alleles(@als); 
    is_deeply(\@ret,["CA","C"],"Testing normalize_alleles(".join(',',@als).")");
    @als = ("CA","CT"); @ret = $vcf->normalize_alleles(@als); 
    is_deeply(\@ret,["CA","CT"],"Testing normalize_alleles(".join(',',@als).")");
    @als = ("GAACCCACA","GA"); @ret = $vcf->normalize_alleles_pos(@als); 
    is_deeply(\@ret,[0,"GAACCCAC","G"],"Testing normalize_alleles_pos(".join(',',@als).")");
    @als = ("CAGTAAAA","CAGAAAA"); @ret = $vcf->normalize_alleles_pos(@als); 
    is_deeply(\@ret,[2,"GT","G"],"Testing normalize_alleles_pos(".join(',',@als).")");
    @als = ("CAGTAAA","CAGAAAA"); @ret = $vcf->normalize_alleles_pos(@als); 
    is_deeply(\@ret,[3,"T","A"],"Testing normalize_alleles_pos(".join(',',@als).")");
    @als = ("GA","GACC"); @ret = $vcf->normalize_alleles_pos(@als); 
    is_deeply(\@ret,[1,"A","ACC"],"Testing normalize_alleles_pos(".join(',',@als).")");
}


