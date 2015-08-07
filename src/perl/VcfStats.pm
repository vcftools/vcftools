#
# Author: petr.danecek@sanger
#

=head1 NAME

VcfStats.pm.  Module for collecting stats from VCF files. 

=head1 SYNOPSIS

    use VcfStats;

    my $vstats = VcfStats->new(file=>'example.vcf.gz');
    while (my $x=$vstats->next_data_hash()) 
    {
        $vstats->collect_stats($x);
    }
    $vstats->dump();

=cut

package VcfStats;
use strict;
use warnings;
use Carp;
use Data::Dumper;
use base 'Vcf';

=head2 new

    About   : Creates new VcfStats.
    Usage   : my $vstats = VcfStats->new(file=>'my.vcf');
    Args    : See Vcf.pm

=cut

sub new
{
    my ($class,@args) = @_;
    my $self = $class->SUPER::new(@args);
    for my $version (@{$$self{versions}})
    {
        if ( $self->isa($version) ) { eval "use base '$version'"; }
    }
    bless($self,$class);
    return $self;
}


sub parse_header
{
    my ($self,@args) = @_;
    $self->SUPER::parse_header(@args);
}

=head2 get_stats_key

    About   : Creates relevant stats hash key, used by select_stats
    Usage   : 
    Args [1]: Hash with filter definition (value to match, range, etc.)
         [2]: Prefix of the stat
         [3]: Value of the filter

=cut

sub get_stats_key
{
    my ($self,$filter,$key,$value) = @_;

    my $stat_key;
    if ( $$filter{exact} )
    {
        if ( $value ne $$filter{value} ) { next; }
        $stat_key = $key.'/'.$value;
    }
    elsif ( $value eq '.' ) { $stat_key = $key.'/.'; }
    elsif ( $$filter{any} ) { $stat_key = $key.'/'.$value; }
    elsif ( $$filter{bin} )
    {
        my $bin = int($value/$$filter{bin_size}) * $$filter{bin_size};
        if ( $bin>$$filter{max} ) { $bin=">$$filter{max}"; }
        $stat_key = $key.'/'.$bin;
    }
    else { $self->throw("TODO: $key...\n"); }
    return $stat_key;
}

=head2 select_stats

    About   : Selects relevant stats hashes
    Usage   : 
    Args [1]: Hash record from next_data_hash
         [2]: Filters

=cut

sub select_stats
{
    my ($self,$rec,$filters) = @_;

    if ( !exists($$self{stats}{all}) ) { $$self{stats}{all}={}; }
    my @mandatory = ( $$self{stats}{all} );

    my %samples;
    for my $sample (keys %{$$rec{gtypes}})
    {
        if ( !exists($$self{stats}{samples}{$sample}) ) 
        { 
            $$self{stats}{samples}{$sample} = {}; 
        }
        push @{$samples{$sample}}, $$self{stats}{samples}{$sample};
    }

    if ( !defined $filters ) { return (\@mandatory,\%samples); }

    while (my ($key,$filter) = each %$filters)
    {
        if ( $key eq 'FILTER' ) 
        { 
            for my $value (@{$$rec{FILTER}})
            {
                my $stats_key = $self->get_stats_key($filter,$key,$value);
                if ( !exists($$self{stats}{$stats_key}) ) { $$self{stats}{$stats_key}={}; }
                push @mandatory, $$self{stats}{$stats_key};
            }
        }
        elsif ( $key eq 'QUAL' ) 
        { 
            my $stats_key = $self->get_stats_key($filter,$key,$$rec{QUAL});
            if ( !exists($$self{stats}{$stats_key}) ) { $$self{stats}{$stats_key}={}; }
            push @mandatory, $$self{stats}{$stats_key};
        }
        elsif ( $key=~m{^INFO/} ) 
        { 
            if ( $$filter{is_flag} )
            {
                if ( $$filter{value} && !exists($$rec{INFO}{$$filter{tag}}) ) { next; }
                elsif ( !$$filter{value} && exists($$rec{INFO}{$$filter{tag}}) ) { next; }
                if ( !exists($$self{stats}{$key}) ) { $$self{stats}{$key}={}; }
                push @mandatory, $$self{stats}{$key};
                next;
            }
            elsif ( exists($$rec{INFO}{$$filter{tag}}) )
            {
                my $stats_key = $self->get_stats_key($filter,$key,$$rec{INFO}{$$filter{tag}});
                if ( !exists($$self{stats}{$stats_key}) ) { $$self{stats}{$stats_key}={}; }
                push @mandatory, $$self{stats}{$stats_key};
            }
        }
        elsif ( $key=~m{^FORMAT/([^/]+)$} )
        {
            while (my ($sample,$hash) = each %{$$rec{gtypes}})
            {
                if ( !exists($$hash{$1}) ) { next; }
                my $stats_key = $self->get_stats_key($filter,$1,$$hash{$1});
                if ( !exists($$self{stats}{samples}{$sample}{user}{$stats_key}) ) { $$self{stats}{samples}{$sample}{user}{$stats_key}={}; }
                push @{$samples{$sample}}, $$self{stats}{samples}{$sample}{user}{$stats_key};
            }
        }
        elsif ( $key=~m{^SAMPLE/([^/]+)/([^/]+)$} )
        {
            if ( !exists($$rec{gtypes}{$1}{$2}) ) { next; }
            my $stats_key = $self->get_stats_key($filter,$2,$$rec{gtypes}{$1}{$2});
            if ( !exists($$self{stats}{samples}{$1}{user}{$stats_key}) ) { $$self{stats}{samples}{$1}{user}{$stats_key}={} }
            push @{$samples{$1}}, $$self{stats}{samples}{$1}{user}{$stats_key};
        }
        else { $self->throw("The feature currently not recognised: $key.\n"); } 
    }
    return (\@mandatory,\%samples);
}

=head2 collect_stats

    About   : Collect stats 
    Usage   : my $x=$vstats->next_data_hash(); $vstats->collect_stats($x);
    Args    : 

=cut

sub collect_stats
{
    my ($self,$rec,$filters) = @_;

    # Ts/Tv and custom numbers based on INFO, QUAL etc. for the mandatory columns
    my ($mandatory_stats,$sample_stats) = $self->select_stats($rec,$filters);
    $self->collect_stats_mandatory($rec,$mandatory_stats);

    # Ts/Tv for samples
    while (my ($sample,$stats) = each %$sample_stats)
    {
        $self->collect_stats_sample($rec,$sample,$stats);
    }

    my %type_keys = ( r=>'ref', s=>'snp', i=>'indel' );

    # Private calls and the number of shared SNPs. Check if:
    #   - there is a nonref variant present only in this sample (samples->sample_name->private)
    #   - there is a nonref variant in N samples (samples->all->shared)
    #   - there is a non-empty call (samples->sample_name->count)
    my $shared = 0;
    my $sample_name;
    while (my ($sample,$stats) = each %$sample_stats)
    {
        my ($alleles,$seps,$is_phased,$is_empty) = $self->parse_haplotype($rec,$sample);
        if ( $is_empty ) { next; }

        my $is_hom=1;
        my %types;
        my $is_ref = 1;
        for my $al (@$alleles)
        {
            if ( $$alleles[0] ne $al ) { $is_hom=0; }
            my ($type,$len,$ht) = $self->event_type($rec,$al);
            $types{$type} = 1;
            if ( $type eq 'r' ) { next; }
            $is_ref = 0;
        }
        for my $stat (@$stats) { $$stat{count}++; }

        for my $type (keys %types)
        {
            my $key = exists($type_keys{$type}) ? $type_keys{$type} : 'other';
            $key .= '_count';
            for my $stat (@$stats) { $$stat{$key}++; }
        }
        my $key;
        if ( exists($types{r}) ) 
        { 
            if ( $is_hom ) { $key='hom_RR'; }
            else { $key='het_RA' }
        }
        elsif ( $is_hom ) { $key='hom_AA'; }
        else { $key='het_AA'; }
        $key .= '_count';
        for my $stat (@$stats) { $$stat{$key}++; }

        $key = $is_phased ? 'phased' : 'unphased';
        for my $stat (@$stats) { $$stat{$key}++; }

        if ( $is_ref ) { next; }
        $shared++;
        if ( !defined $sample_name ) { $sample_name = $sample; }
    }
    $$self{stats}{all}{shared}{$shared}++;
    if ( $shared==1 )
    {
        for my $stat (@{$$sample_stats{$sample_name}}) { $$stat{private}++; }
    }
}


=head2 collect_stats_mandatory

    About   : Collect stats based on mandatory columns
    Usage   : my $x=$vstats->next_data_hash(); $vstats->collect_stats_mandatory($x);
    Args    : 

=cut

sub collect_stats_mandatory
{
    my ($self,$rec,$stats) = @_;

    # How many mono,bi,tri-allelic etc sites are there
    my $nalt = 0;
    if ( !scalar keys %{$$rec{gtypes}} ) 
    { 
        $nalt = scalar @{$$rec{ALT}};
        if ( $nalt==1 && $$rec{ALT}[0] eq '.' ) { $nalt=0 }
    }
    elsif ( exists($$rec{INFO}{AC}) )
    {
        for my $ac (split(/,/,$$rec{INFO}{AC})) 
        { 
            if ( $ac ) { $nalt++; } 
        }
    }
    else
    {
        my ($an,$ac,$acs) = $self->calc_an_ac($$rec{gtypes});
        for my $ac (@$acs)
        {
            if ( $ac ) { $nalt++; }
        }
    }

    my %types;
    for my $alt (@{$$rec{ALT}})
    {
        if ( $alt eq '.' ) { $alt=$$rec{REF}; }
        my $type = $self->add_variant($rec,$alt,$stats);
        $types{$type} = 1;
    }


    # Increment counters
    for my $stat (@$stats)
    {
        $$stat{'nalt_'.$nalt}++;
        $$stat{count}++;
        for my $type (keys %types)
        {
            $$stat{$type.'_count'}++;
        }
    }
}


=head2 collect_stats_sample

    About   : Collect stats for given sample
    Usage   : my $x=$vstats->next_data_hash(); $vstats->collect_stats_sample($x,'NA0001');
    Args [1]  hash row from next_data_hash 
         [2]  sample name
         [3]  stats to collect

=cut

sub collect_stats_sample
{
    my ($self,$rec,$sample,$stats) = @_;

    my ($alleles,$seps,$is_phased,$is_empty) = $self->parse_haplotype($rec,$sample);
    if ( @$alleles > 2 ) { $self->throw("FIXME: currently handling diploid data only (easy to fix)\n"); }
    my $prev;
    for my $al (@$alleles)
    {
        if ( !defined $prev or $prev ne $al )
        {
            # Only heterozygous SNPs will be counted twice
            $self->add_variant($rec,$al,$stats);
        }
        $prev = $al;
    }
}


=head2 add_variant

    About   : Register mutation type in the selected pool
    Usage   : $vstats->add_variant('A','AT',$stats);
              $vstats->add_variant($rec,'AT',$stats);
    Args      [1] Reference haplotype or VCF data line parsed by next_data_hash
              [2] Variant haplotype
              [3] Array of hash stats
    Returns : The event type (snp,indel,ref)

=cut

sub add_variant
{
    my ($self,$ref,$alt,$stats) = @_;
    my $key_type = 'other';
    my %key_subt;
    if ( $alt eq '.' )
    {
        $key_type = 'missing';
    }
    else
    {
        my ($type,$len,$ht) = $self->event_type($ref,$alt);
        if ( $type eq 's' ) 
        { 
            $key_type = 'snp';

            # The SNP can be encoded for example as GTTTTTTT>CTTTTTTT
            my $ref_str = ref($ref) eq 'HASH' ? $$ref{REF} : $ref;
            my $ref_len = length($ref_str);
            if ( $ref_len>1 )
            {
                for (my $i=0; $i<$ref_len; $i++)
                {
                    my $ref_nt = substr($ref_str,$i,1);
                    my $alt_nt = substr($alt,$i,1);
                    if ( $ref_nt ne $alt_nt )
                    {
                        $key_subt{$ref_nt.'>'.$alt_nt}++;
                    }
                }
            }
            else
            {
                $key_subt{$ref_str.'>'.$alt}++;
            }
        }
        elsif ( $type eq 'i' ) { $key_type = 'indel'; $key_subt{$len}++; }
        elsif ( $type eq 'r' ) { $key_type = 'ref'; }
    }
    for my $stat (@$stats)
    {
        if ( %key_subt )
        {
            while (my ($subt,$value)=each %key_subt)
            {
                $$stat{$key_type}{$subt}+=$value;
            }
        }
        else
        {
            $$stat{$key_type}++;
        }
    }
    return $key_type;
}


=head2 dump

    About   : Produce Data::Dumper dump of the collected stats
    Usage   : 
    Args    :
    Returns : The dump.

=cut

sub dump
{
    my ($self) = @_;
    return Dumper($$self{stats});
}


sub _calc_tstv
{
    my ($self,$stat) = @_;

    my $ts = 0;
    for my $mut (qw(A>G G>A C>T T>C))
    {
        if ( exists($$stat{$mut}) ) { $ts += $$stat{$mut}; }
    }
    my $tv = 0;
    for my $mut (qw(A>C C>A G>T T>G A>T T>A C>G G>C))
    {
        if ( exists($$stat{$mut}) ) { $tv += $$stat{$mut}; }
    }
    my $ratio = $tv ? $ts/$tv : 0;
    return ($ts,$tv,$ratio);
}


=head2 dump_tstv

    About   : Calculate transitions/transversions ratio and output string
    Usage   : 
    Args    :
    Returns : Formatted string

=cut

sub dump_tstv
{
    my ($self,$stats) = @_;
    my $out = "#Transitions\tTransversions\tts/tv\tSample\n";
    for my $key (sort keys %$stats)
    {
        if ( !exists($$stats{$key}{snp}) ) { next; }
        my $stat = $$stats{$key}{snp};
        my ($ts,$tv,$ratio) = $self->_calc_tstv($stat);
        $out .= sprintf "%d\t%d\t%.2f\t%s\n", $ts,$tv,$ratio,$key;
    }
    return $out;
}


=head2 dump_qual_tstv

    About   : Calculate marginal transitions/transversions ratios for QUAL/* stats
    Usage   : 
    Args    :
    Returns : Formatted string

=cut

sub dump_qual_tstv
{
    my ($self,$file) = @_;
    my @values;
    for my $stat (keys %{$$self{stats}})
    {
        if ( !($stat=~m{^QUAL/(.+)}) ) { next; }
        my $qual  = $1;
        # The quality record can be also of the form ">200". Exclude these from numeric comparison
        if ( !($qual=~/^[0-9.]+$/) ) { $qual = "#$qual"; }
        my $count = $$self{stats}{$stat}{count};
        if ( !exists($$self{stats}{$stat}{snp}) ) { next; }
        my ($ts,$tv,$ratio) = $self->_calc_tstv($$self{stats}{$stat}{snp});
        push @values, [$qual,$count,$ratio];
    }
    my @svalues = sort { if ($$a[0]=~/^#/ or $$b[0]=~/^#/) { return $$a[0] cmp $$b[0]; } return $$a[0] <=> $$b[0]; } @values;
    my $out = "#Quality\tMarginal count\tMarginal Ts/Tv\n";
    for my $val (@svalues)
    {
        if ( $$val[0]=~/^#/ )
        {
            $out .= sprintf "%s\t%d\t%.2f\n", $$val[0],$$val[1],$$val[2];
        }
        else
        {
            $out .= sprintf "%.2f\t%d\t%.2f\n", $$val[0],$$val[1],$$val[2];
        }
    }
    return $out;
}


=head2 dump_counts

    About   : 
    Usage   : 
    Args    :
    Returns : Formatted string

=cut

sub dump_counts
{
    my ($self) = @_;
    my $out = "#Count\tFilter\n";
    for my $key (sort keys %{$$self{stats}})
    {
        if ( !exists($$self{stats}{$key}{count}) ) { next; }
        $out .= sprintf "%d\t%s\n", $$self{stats}{$key}{count},$key;
    }
    for my $key (sort keys %{$$self{stats}{samples}})
    {
        if ( !exists($$self{stats}{samples}{$key}{count}) ) { next; }
        $out .= sprintf "%d\tsamples/%s\n", $$self{stats}{samples}{$key}{count},$key;
    }
    return $out;
}

sub dump_snp_counts
{
    my ($self) = @_;
    my $out = "#Count\tFilter\n";
    for my $key (sort keys %{$$self{stats}})
    {
        if ( !exists($$self{stats}{$key}{snp_count}) ) { next; }
        $out .= sprintf "%d\t%s\n",$$self{stats}{$key}{snp_count},$key;
    }
    for my $key (sort keys %{$$self{stats}{samples}})
    {
        if ( !exists($$self{stats}{samples}{$key}{snp_count}) ) { next; }
        $out .= sprintf "%d\tsamples/%s\n", $$self{stats}{samples}{$key}{snp_count},$key;
    }
    return $out;
}

sub dump_indel_counts
{
    my ($self) = @_;
    my $out = "#Count\tFilter\n";
    for my $key (sort keys %{$$self{stats}})
    {
        if ( !exists($$self{stats}{$key}{indel_count}) ) { next; }
        $out .= sprintf "%d\t%s\n",$$self{stats}{$key}{indel_count},$key;
    }
    for my $key (sort keys %{$$self{stats}{samples}})
    {
        if ( !exists($$self{stats}{samples}{$key}{indel_count}) ) { next; }
        $out .= sprintf "%d\tsamples/%s\n", $$self{stats}{samples}{$key}{indel_count},$key;
    }
    return $out;
}

sub dump_shared_counts
{
    my ($self) = @_;
    my $out = "#Shared SNPs\tFrequency\n";
    for my $key (sort {$a<=>$b} keys %{$$self{stats}{all}{shared}})
    {
        $out .= sprintf "%d\t%s\n", $key,$$self{stats}{all}{shared}{$key};
    }
    return $out;
}

sub dump_private_counts
{
    my ($self) = @_;
    my $out = "#Private SNPs\tSample\n";
    for my $key (sort keys %{$$self{stats}{samples}})
    {
        if ( !exists($$self{stats}{samples}{$key}{private}) ) { next; }
        $out .= sprintf "%d\t%s\n", $$self{stats}{samples}{$key}{private},$key;
    }
    return $out;
}

sub _init_path
{
    my ($self,$prefix) = @_;
    if ( $prefix=~m{/} )
    {
        # A directory should be created. This will populate dir and prefix, for example
        #   prefix  -> dir      prefix
        #   ----------------------------
        #   out                 out.dump
        #   out/       out/     out/out.dump
        #   out/xxx    out/     out/xxx.dump 
        #
        my $dir = '';
        if ( $prefix=~m{/[^/]+$} ) { $dir=$`; }
        elsif ( $prefix=~m{/([^/]+)/$} ) { $dir = $`.'/'.$1; $prefix = $dir.'/'.$1; }
        elsif ( $prefix=~m{([^/]+)/?$} ) { $dir=$1; $prefix=$dir.'/'.$1; }
        if ( $dir ) { `mkdir -p $dir`; }
    }
    return $prefix;
}

sub legend
{
    my ($self) = @_;
    return q[
count
    Number of positions with known genotype

nalt_X
    Number of monoallelic (X=0), biallelic (X=1), etc. sites
    
ref, ref_count
    Number of sites containing reference allele

shared
    Number of sites having a non-reference allele in 0,1,2,etc samples

snp_count
    Number of positions with SNPs
    ];
}



=head2 save_stats

    About   : Save all collected stats to files
    Usage   : 
    Args    : The prefix of output files. Non-existent directories will be created.
    Returns : N/A

=cut

sub save_stats
{
    my ($self,$prefix) = @_;

    if ( !defined $prefix )
    {
        print $self->dump();
        return;
    }

    my $path = $self->_init_path($prefix);
    $self->_write_file($path.'.legend', $self->legend());
    $self->_write_file($path.'.dump', $self->dump());
    $self->_write_file($path.'.tstv', $self->dump_tstv($$self{stats}));
    $self->_write_file($path.'.counts', $self->dump_counts());
    $self->_write_file($path.'.snps', $self->dump_snp_counts());
    $self->_write_file($path.'.indels', $self->dump_indel_counts());
    $self->_write_file($path.'.qual-tstv',$self->dump_qual_tstv);
    $self->_write_file($path.'.shared',$self->dump_shared_counts());
    $self->_write_file($path.'.private',$self->dump_private_counts());

    if ( exists($$self{stats}{samples}) )
    {
        $self->_write_file($path.'.samples-tstv',$self->dump_tstv($$self{stats}{samples}));
    }
}


sub _write_file
{
    my ($self,$fname,$text) = @_;
    open(my $fh,'>',$fname) or $self->throw("$fname: $!");
    print $fh $text;
    close($fh);
}


1;

