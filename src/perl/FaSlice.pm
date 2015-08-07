# Author: petr.danecek@sanger
#

=head1 NAME

FaSlice.pm.  Module for cached access to fasta sequences, employs samtools faidx. 

=head1 SYNOPSIS

use FaSlice;

my $fa = FaSlice->new(file=>'ref.fa');
$fa->get_base(1,12345);
$fa->get_slice(1,12345,54321);

=cut


package FaSlice;

use strict;
use warnings;
use Carp;

=head2 new

    About   : Creates new FaSlice object.
    Usage   : my $fa = FaSlice->new(file=>'ref.fa');
    Args    : file   .. the fasta file
              oob    .. out-of-bounds requests: one of 'throw' (throws), 'N' (fills the missing bases with Ns), or '' (returns empty string, default)
              size   .. size of the cached chunk read by samtools faidx (1_000_000)

=cut

sub new
{
    my ($class,@args) = @_;
    my $self = @args ? {@args} : {};
    bless $self, ref($class) || $class;
    if ( !$$self{file} ) { $self->throw("Missing the parameter file\n"); }
    $$self{chr}  = undef;
    $$self{from} = undef;
    $$self{to}   = undef;
    if ( !$$self{size} ) { $$self{size}=1_000_000; }
    $$self{ncache_missed} = 0;
    $$self{nqueries} = 0;
    if ( !exists($$self{oob}) ) { $$self{oob}=''; }
    if ( $$self{oob} ne '' && $$self{oob} ne 'throw' && $$self{oob} ne 'N' ) { $self->throw("The value of oob not recognised: [$$self{oob}]"); }
    $self->chromosome_naming($$self{file});
    return $self;
}

sub throw
{
    my ($self,@msg) = @_;
    confess(@msg);
}

sub cmd
{
    my ($self,$cmd) = @_;
    my @out = `$cmd`;
    if ( $? )
    {
        my @msg = ();
        push @msg, qq[The command "$cmd" returned non-zero status $?];
        if ( $! )
        {
            push @msg, ": $!\n";
        }
        else
        {
            push @msg, ".\n";
        }
        if ( scalar @out )
        {
            push @msg, @out;
        }
        $self->throw(@msg);
    }
    return (@out);
}

# Read the first file of the fasta file and make a guess: Are all chromosomes
#   names as 'chr1','chr2',etc or just '1','2',...?
# Future TODO: more robust chromosome name mapping?
sub chromosome_naming
{
    my ($self,$fa_file) = @_;
    open(my $fh,'<',"$fa_file.fai") or $self->throw("$fa_file.fai: $!");
    my $line=<$fh>;
    if ( !($line=~/^(chr)?\S+\t/) ) { chomp($line); $self->throw("FIXME: the sequence names not in '>(chr)?\\S+' format [$line] ... $fa_file.fai\n"); }
    close($fh); 
    $$self{chr_naming} = defined $1 ? $1 : '';
}

sub cache_chr_lengths
{
    my ($self) = @_;
    if ( exists($$self{chr_lengths}) ) { return; }
    open(my $fh,'<',"$$self{file}.fai") or $self->throw("$$self{file}.fai: $!");
    while (my $line=<$fh>)
    {
        my @items = split(/\t/,$line);
        my $chr = $$self{chr_naming}.$items[0];
        $$self{chr_lengths}{$chr} = $items[1];
    }
    close($fh) or $self->throw("close $$self{file}.fai");
}

sub read_chunk
{
    my ($self,$chr,$pos) = @_;
    $$self{chr}  = $chr;
    $chr =~ s/^chr//;
    $chr = $$self{chr_naming}.$chr;
    if ( exists($$self{chr_lengths}) && (!exists($$self{chr_lengths}{$chr}) or $$self{chr_lengths}{$chr} < $pos ) )
    {
        $$self{to} = $$self{from} - 1;
        $$self{chunk} = '';
        return;
    }
    my $to = $pos + $$self{size};
    my $cmd = "samtools faidx $$self{file} $chr:$pos-$to";
    my @out = $self->cmd($cmd) or $self->throw("$cmd: $!");
    my $line = shift(@out);
    if ( !($line=~/^>$chr:(\d+)-(\d+)/) ) { $self->throw("Could not parse: $line"); }
    $$self{from} = $1;
    my $chunk = '';
    while ($line=shift(@out))
    {
        chomp($line);
        $chunk .= $line;
    }
    $$self{to} = $$self{from} + length($chunk) - 1;
    $$self{chunk} = $chunk;
    $self->cache_chr_lengths();
    return;
}

=head2 get_base

    About   : Retrieves base at the given chromosome and position
    Usage   : my $fa = FaSlice->new(file=>'ref.fa'); $fa->get_base(1,12345);
    Args    : chromosome
              1-based coordinate

=cut

sub get_base
{
    my ($self,$chr,$pos) = @_;
    if ( !$$self{chr} || $chr ne $$self{chr} || $pos<$$self{from} || $pos>$$self{to} )
    {
        $self->read_chunk($chr,$pos);
    }
    $$self{nqueries}++;
    my $idx = $pos - $$self{from};
    if ( $$self{from}>$$self{to} ) 
    { 
        if ( $$self{oob} eq '' ) { return ''; }
        elsif ( $$self{oob} eq 'N' ) { return 'N'; }
        $self->throw("No such site $chr:$pos in $$self{file}\n"); 
    }
    return substr($$self{chunk},$idx,1);
}


=head2 get_slice

    About   : Retrieves region
    Usage   : my $fa = FaSlice->new(file=>'ref.fa'); $fa->get_base(1,12345,54321);
    Args    : chromosome
              1-based coordinate

=cut

sub get_slice
{
    my ($self,$chr,$from,$to) = @_;
    if ( $to-$from >= $$self{size} ) { $$self{size} = $to-$from+1; }
    if ( $from>$to ) { $self->throw("Expected $from>$to\n"); }
    if ( !$$self{chr} || $chr ne $$self{chr} || $from<$$self{from} || $to>$$self{to} )
    {
        $self->read_chunk($chr,$from);
    }
    $$self{nqueries}++;

    if ( $$self{from}>$$self{to} || $$self{from}>$from || $$self{to}<$to )
    {
        if ( $$self{oob} eq 'throw' ) { $self->throw("The region out of bounds $chr:$from-$to in $$self{file}\n"); }
        elsif ( $$self{oob} eq '' ) { return ''; }

        if ( $$self{from}>$$self{to} ) { return 'N' x ($to-$from+1); }
        if ( $$self{from}>$to ) { $self->throw("FIXME: this shouldn't happen $chr:$from-$to .. $$self{from},$$self{to} .. $$self{file}"); }

        my $lfill = '';
        my $rfill = '';
        if ( $$self{from}>$from ) { $lfill = 'N' x ($$self{from}-$from); $from=$$self{from}; }
        if ( $$self{to}<$to ) { $rfill = 'N' x ($to-$$self{to}); $to=$$self{to}; }
        return $lfill . substr($$self{chunk},$from-$$self{from},$to-$from+1) . $rfill;
    }
    return substr($$self{chunk},$from-$$self{from},$to-$from+1);
}


# http://www.illumina.com/documents/products/technotes/technote_topbot.pdf
sub illumina_alleles_TOP_to_ref
{
    my ($self,$a1,$a2,$chr,$pos,$ref) = @_;
    my %map = (A=>'T', C=>'G', G=>'C', T=>'A');
    my %top = ( 
            A=>{A=>-2,C=> 1,G=> 1,T=>-1}, 
            C=>{A=> 1,C=>-2,G=>-1,T=> 0}, 
            G=>{A=> 1,C=>-1,G=>-2,T=> 0}, 
            T=>{A=>-1,C=> 0,G=> 0,T=>-2} ); 

    my $stat = $top{$a1}{$a2};
    if ( $stat==-2 ) { $self->throw("Expected two different bases, got $a1 and $a2.\n"); }
    if ( $stat==-1 )
    {
        # Now we should do the sequence walking to see if the reference is TOP or BOT,
        #   but we do not this in ill-to-vcf: C/G would become G/C and A/T would become T/A.
        return ($a1,$a2);
    }
    if ( $stat==0 ) { $self->throw("Expected Illumina TOP, got $a1 and $a2.\n"); }
    if ( $ref eq $a1 or $ref eq $a2 ) { return ($a1,$a2); }
    return ($map{$a1},$map{$a2});
}


1;

