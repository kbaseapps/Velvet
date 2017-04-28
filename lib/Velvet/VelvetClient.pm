package Velvet::VelvetClient;

use JSON::RPC::Client;
use POSIX;
use strict;
use Data::Dumper;
use URI;
use Bio::KBase::Exceptions;
my $get_time = sub { time, 0 };
eval {
    require Time::HiRes;
    $get_time = sub { Time::HiRes::gettimeofday() };
};

use Bio::KBase::AuthToken;

# Client version should match Impl version
# This is a Semantic Version number,
# http://semver.org
our $VERSION = "0.1.0";

=head1 NAME

Velvet::VelvetClient

=head1 DESCRIPTION


A KBase module: Velvet

This is a KBase module that wraps the open source package "Short read de novo assembler using de Bruijn graphs"
Version 1.2.10

References:
https://github.com/dzerbino/velvet
https://github.com/dzerbino/velvet/blob/master/Columbus_manual.pdf


=cut

sub new
{
    my($class, $url, @args) = @_;
    

    my $self = {
	client => Velvet::VelvetClient::RpcClient->new,
	url => $url,
	headers => [],
    };

    chomp($self->{hostname} = `hostname`);
    $self->{hostname} ||= 'unknown-host';

    #
    # Set up for propagating KBRPC_TAG and KBRPC_METADATA environment variables through
    # to invoked services. If these values are not set, we create a new tag
    # and a metadata field with basic information about the invoking script.
    #
    if ($ENV{KBRPC_TAG})
    {
	$self->{kbrpc_tag} = $ENV{KBRPC_TAG};
    }
    else
    {
	my ($t, $us) = &$get_time();
	$us = sprintf("%06d", $us);
	my $ts = strftime("%Y-%m-%dT%H:%M:%S.${us}Z", gmtime $t);
	$self->{kbrpc_tag} = "C:$0:$self->{hostname}:$$:$ts";
    }
    push(@{$self->{headers}}, 'Kbrpc-Tag', $self->{kbrpc_tag});

    if ($ENV{KBRPC_METADATA})
    {
	$self->{kbrpc_metadata} = $ENV{KBRPC_METADATA};
	push(@{$self->{headers}}, 'Kbrpc-Metadata', $self->{kbrpc_metadata});
    }

    if ($ENV{KBRPC_ERROR_DEST})
    {
	$self->{kbrpc_error_dest} = $ENV{KBRPC_ERROR_DEST};
	push(@{$self->{headers}}, 'Kbrpc-Errordest', $self->{kbrpc_error_dest});
    }

    #
    # This module requires authentication.
    #
    # We create an auth token, passing through the arguments that we were (hopefully) given.

    {
	my %arg_hash2 = @args;
	if (exists $arg_hash2{"token"}) {
	    $self->{token} = $arg_hash2{"token"};
	} elsif (exists $arg_hash2{"user_id"}) {
	    my $token = Bio::KBase::AuthToken->new(@args);
	    if (!$token->error_message) {
	        $self->{token} = $token->token;
	    }
	}
	
	if (exists $self->{token})
	{
	    $self->{client}->{token} = $self->{token};
	}
    }

    my $ua = $self->{client}->ua;	 
    my $timeout = $ENV{CDMI_TIMEOUT} || (30 * 60);	 
    $ua->timeout($timeout);
    bless $self, $class;
    #    $self->_validate_version();
    return $self;
}




=head2 run_velveth

  $output = $obj->run_velveth($params)

=over 4

=item Parameter and return types

=begin html

<pre>
$params is a Velvet.VelvethParams
$output is a Velvet.VelvethResults
VelvethParams is a reference to a hash where the following keys are defined:
	out_folder has a value which is a string
	hash_length has a value which is an int
	filename has a value which is a string
	file_format has a value which is a string
	read_type has a value which is a string
VelvethResults is a reference to a hash where the following keys are defined:
	out_folder has a value which is a string

</pre>

=end html

=begin text

$params is a Velvet.VelvethParams
$output is a Velvet.VelvethResults
VelvethParams is a reference to a hash where the following keys are defined:
	out_folder has a value which is a string
	hash_length has a value which is an int
	filename has a value which is a string
	file_format has a value which is a string
	read_type has a value which is a string
VelvethResults is a reference to a hash where the following keys are defined:
	out_folder has a value which is a string


=end text

=item Description

Definition of run_velveth

=back

=cut

 sub run_velveth
{
    my($self, @args) = @_;

# Authentication: required

    if ((my $n = @args) != 1)
    {
	Bio::KBase::Exceptions::ArgumentValidationError->throw(error =>
							       "Invalid argument count for function run_velveth (received $n, expecting 1)");
    }
    {
	my($params) = @args;

	my @_bad_arguments;
        (ref($params) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument 1 \"params\" (value was \"$params\")");
        if (@_bad_arguments) {
	    my $msg = "Invalid arguments passed to run_velveth:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	    Bio::KBase::Exceptions::ArgumentValidationError->throw(error => $msg,
								   method_name => 'run_velveth');
	}
    }

    my $url = $self->{url};
    my $result = $self->{client}->call($url, $self->{headers}, {
	    method => "Velvet.run_velveth",
	    params => \@args,
    });
    if ($result) {
	if ($result->is_error) {
	    Bio::KBase::Exceptions::JSONRPC->throw(error => $result->error_message,
					       code => $result->content->{error}->{code},
					       method_name => 'run_velveth',
					       data => $result->content->{error}->{error} # JSON::RPC::ReturnObject only supports JSONRPC 1.1 or 1.O
					      );
	} else {
	    return wantarray ? @{$result->result} : $result->result->[0];
	}
    } else {
        Bio::KBase::Exceptions::HTTP->throw(error => "Error invoking method run_velveth",
					    status_line => $self->{client}->status_line,
					    method_name => 'run_velveth',
				       );
    }
}
 


=head2 run_velvetg

  $output = $obj->run_velvetg($params)

=over 4

=item Parameter and return types

=begin html

<pre>
$params is a Velvet.VelvetgParams
$output is a Velvet.VelvetgResults
VelvetgParams is a reference to a hash where the following keys are defined:
	wk_folder has a value which is a string
	cov_cutoff has a value which is a float
	ins_length has a value which is an int
	read_trkg has a value which is an int
	min_contig_length has a value which is an int
	amos_file has a value which is an int
	exp_cov has a value which is a float
	long_cov_cutoff has a value which is a float
VelvetgResults is a reference to a hash where the following keys are defined:
	wk_folder has a value which is a string

</pre>

=end html

=begin text

$params is a Velvet.VelvetgParams
$output is a Velvet.VelvetgResults
VelvetgParams is a reference to a hash where the following keys are defined:
	wk_folder has a value which is a string
	cov_cutoff has a value which is a float
	ins_length has a value which is an int
	read_trkg has a value which is an int
	min_contig_length has a value which is an int
	amos_file has a value which is an int
	exp_cov has a value which is a float
	long_cov_cutoff has a value which is a float
VelvetgResults is a reference to a hash where the following keys are defined:
	wk_folder has a value which is a string


=end text

=item Description

Definition of run_velvetg

=back

=cut

 sub run_velvetg
{
    my($self, @args) = @_;

# Authentication: required

    if ((my $n = @args) != 1)
    {
	Bio::KBase::Exceptions::ArgumentValidationError->throw(error =>
							       "Invalid argument count for function run_velvetg (received $n, expecting 1)");
    }
    {
	my($params) = @args;

	my @_bad_arguments;
        (ref($params) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument 1 \"params\" (value was \"$params\")");
        if (@_bad_arguments) {
	    my $msg = "Invalid arguments passed to run_velvetg:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	    Bio::KBase::Exceptions::ArgumentValidationError->throw(error => $msg,
								   method_name => 'run_velvetg');
	}
    }

    my $url = $self->{url};
    my $result = $self->{client}->call($url, $self->{headers}, {
	    method => "Velvet.run_velvetg",
	    params => \@args,
    });
    if ($result) {
	if ($result->is_error) {
	    Bio::KBase::Exceptions::JSONRPC->throw(error => $result->error_message,
					       code => $result->content->{error}->{code},
					       method_name => 'run_velvetg',
					       data => $result->content->{error}->{error} # JSON::RPC::ReturnObject only supports JSONRPC 1.1 or 1.O
					      );
	} else {
	    return wantarray ? @{$result->result} : $result->result->[0];
	}
    } else {
        Bio::KBase::Exceptions::HTTP->throw(error => "Error invoking method run_velvetg",
					    status_line => $self->{client}->status_line,
					    method_name => 'run_velvetg',
				       );
    }
}
 
  
sub status
{
    my($self, @args) = @_;
    if ((my $n = @args) != 0) {
        Bio::KBase::Exceptions::ArgumentValidationError->throw(error =>
                                   "Invalid argument count for function status (received $n, expecting 0)");
    }
    my $url = $self->{url};
    my $result = $self->{client}->call($url, $self->{headers}, {
        method => "Velvet.status",
        params => \@args,
    });
    if ($result) {
        if ($result->is_error) {
            Bio::KBase::Exceptions::JSONRPC->throw(error => $result->error_message,
                           code => $result->content->{error}->{code},
                           method_name => 'status',
                           data => $result->content->{error}->{error} # JSON::RPC::ReturnObject only supports JSONRPC 1.1 or 1.O
                          );
        } else {
            return wantarray ? @{$result->result} : $result->result->[0];
        }
    } else {
        Bio::KBase::Exceptions::HTTP->throw(error => "Error invoking method status",
                        status_line => $self->{client}->status_line,
                        method_name => 'status',
                       );
    }
}
   

sub version {
    my ($self) = @_;
    my $result = $self->{client}->call($self->{url}, $self->{headers}, {
        method => "Velvet.version",
        params => [],
    });
    if ($result) {
        if ($result->is_error) {
            Bio::KBase::Exceptions::JSONRPC->throw(
                error => $result->error_message,
                code => $result->content->{code},
                method_name => 'run_velvetg',
            );
        } else {
            return wantarray ? @{$result->result} : $result->result->[0];
        }
    } else {
        Bio::KBase::Exceptions::HTTP->throw(
            error => "Error invoking method run_velvetg",
            status_line => $self->{client}->status_line,
            method_name => 'run_velvetg',
        );
    }
}

sub _validate_version {
    my ($self) = @_;
    my $svr_version = $self->version();
    my $client_version = $VERSION;
    my ($cMajor, $cMinor) = split(/\./, $client_version);
    my ($sMajor, $sMinor) = split(/\./, $svr_version);
    if ($sMajor != $cMajor) {
        Bio::KBase::Exceptions::ClientServerIncompatible->throw(
            error => "Major version numbers differ.",
            server_version => $svr_version,
            client_version => $client_version
        );
    }
    if ($sMinor < $cMinor) {
        Bio::KBase::Exceptions::ClientServerIncompatible->throw(
            error => "Client minor version greater than Server minor version.",
            server_version => $svr_version,
            client_version => $client_version
        );
    }
    if ($sMinor > $cMinor) {
        warn "New client version available for Velvet::VelvetClient\n";
    }
    if ($sMajor == 0) {
        warn "Velvet::VelvetClient version is $svr_version. API subject to change.\n";
    }
}

=head1 TYPES



=head2 report_ref

=over 4



=item Description

A 'typedef' allows you to provide a more specific name for
a type.  Built-in primitive types include 'string', 'int',
'float'.  Here we define a type named assembly_ref to indicate
a string that should be set to a KBase ID reference to an
Assembly data object.


=item Definition

=begin html

<pre>
a string
</pre>

=end html

=begin text

a string

=end text

=back



=head2 VelvethParams

=over 4



=item Description

Arguments for run_velveth
velveth, help

Compilation settings:
CATEGORIES = 2
MAXKMERLENGTH = 31

Usage:
./velveth directory hash_length {[-file_format][-read_type][-separate|-interleaved] filename1 [filename2 ...]} {...} [options]

        directory       : directory name for output files
        hash_length     : EITHER an odd integer (if even, it will be decremented) <= 31 (if above, will be reduced)
                        : OR: m,M,s where m and M are odd integers (if not, they will be decremented) with m < M <= 31 (if above, will be reduced)
                                and s is a step (even number). Velvet will then hash from k=m to k=M with a step of s
        filename        : path to sequence file or - for standard input

File format options:
        -fasta  -fastq  -raw    -fasta.gz       -fastq.gz       -raw.gz -sam    -bam    -fmtAuto
        (Note: -fmtAuto will detect fasta or fastq, and will try the following programs for decompression : gunzip, pbunzip2, bunzip2

File layout options for paired reads (only for fasta and fastq formats):
        -interleaved    : File contains paired reads interleaved in the one file (default)
        -separate       : Read 2 separate files for paired reads

Read type options:
        -short  -shortPaired
        -short2 -shortPaired2
        -long   -longPaired
        -reference

Options:
        -strand_specific        : for strand specific transcriptome sequencing data (default: off)
        -reuse_Sequences        : reuse Sequences file (or link) already in directory (no need to provide original filenames in this case (default: off)
        -reuse_binary   : reuse binary sequences file (or link) already in directory (no need to provide original filenames in this case (default: off)
        -noHash                 : simply prepare Sequences file, do not hash reads or prepare Roadmaps file (default: off)
        -create_binary          : create binary CnyUnifiedSeq file (default: off)

Synopsis:

- Short single end reads:
        velveth Assem 29 -short -fastq s_1_sequence.txt

- Paired-end short reads (remember to interleave paired reads):
        velveth Assem 31 -shortPaired -fasta interleaved.fna

- Paired-end short reads (using separate files for the paired reads)
        velveth Assem 31 -shortPaired -fasta -separate left.fa right.fa

- Two channels and some long reads:
        velveth Assem 43 -short -fastq unmapped.fna -longPaired -fasta SangerReads.fasta

- Three channels:
        velveth Assem 35 -shortPaired -fasta pe_lib1.fasta -shortPaired2 pe_lib2.fasta -short3 se_lib1.fa

Output:
        out_folder/Roadmaps
        out_folder/Sequences
                [Both files are picked up by graph, so please leave them there]

    Here is the test examples and their stdout printouts:
root@c50eaaa56231:/kb/module# ls /data
__READY__  velvet_data
root@c50eaaa56231:/kb/module# cd /velvet_data/
root@c50eaaa56231:/velvet_data# ls
test_long.fa  test_reads.fa  test_reads.sam  test_reference.fa
root@c50eaaa56231:/velvet_data# sort test_reads.sam > mySortedReads.sam
root@c50eaaa56231:/velvet_data# velveth test_dir 21 -reference test_reference.fa -shortPaired -sam mySortedReads.sam
[0.000000] Reading FastA file test_reference.fa;
[0.006270] 1 sequences found
[0.006299] Done
[0.006331] Reading SAM file mySortedReads.sam
[0.246146] 142858 reads found.
[0.246170] Done
[0.246172] Reference mapping counters
[0.246173] Name Read mappings
[0.246174] SEQUENCE     142858
[0.246222] Reading read set file test_dir/Sequences;
[0.259455] 142859 sequences found
[0.259559] Read 1 of length 32773, longer than limit 32767
[0.259575] You should modify recompile with the LONGSEQUENCES option (cf. manual)
        string out_folder; #folder name for output files
        int hash_length; #EITHER an odd integer (if even, it will be decremented) <= 31 (if above, will be reduced)L
        string filename; #path to sequence file or - for standard input
        string file_format; #e.g., -fasta, -fastq, -raw,-fasta.gz, -fastq.gz, -raw.gz, -sam, -bam, -fmtAuto
        string read_type; #e.g., -short (-shortPaired), -long(-longPaired), or -reference


=item Definition

=begin html

<pre>
a reference to a hash where the following keys are defined:
out_folder has a value which is a string
hash_length has a value which is an int
filename has a value which is a string
file_format has a value which is a string
read_type has a value which is a string

</pre>

=end html

=begin text

a reference to a hash where the following keys are defined:
out_folder has a value which is a string
hash_length has a value which is an int
filename has a value which is a string
file_format has a value which is a string
read_type has a value which is a string


=end text

=back



=head2 VelvethResults

=over 4



=item Definition

=begin html

<pre>
a reference to a hash where the following keys are defined:
out_folder has a value which is a string

</pre>

=end html

=begin text

a reference to a hash where the following keys are defined:
out_folder has a value which is a string


=end text

=back



=head2 VelvetgParams

=over 4



=item Description

Arguments for run_velvetg
wk_folder                       : working directory name

Standard options:
-cov_cutoff <floating-point|auto>       : removal of low coverage nodes AFTER tour bus or allow the system to infer it
        (default: no removal)
-ins_length <integer>           : expected distance between two paired end reads (default: no read pairing)
-read_trkg <yes|no>             : tracking of short read positions in assembly (default: no tracking)
-min_contig_lgth <integer>      : minimum contig length exported to contigs.fa file (default: hash length * 2)
-amos_file <yes|no>             : export assembly to AMOS file (default: no export)
-exp_cov <floating point|auto>  : expected coverage of unique regions or allow the system to infer it
        (default: no long or paired-end read resolution)
-long_cov_cutoff <floating-point>: removal of nodes with low long-read coverage AFTER tour bus
        (default: no removal)

Advanced options:
-ins_length* <integer>          : expected distance between two paired-end reads in the respective short-read dataset (default: no read pairing)
-ins_length_long <integer>      : expected distance between two long paired-end reads (default: no read pairing)
-ins_length*_sd <integer>       : est. standard deviation of respective dataset (default: 10% of corresponding length)
        [replace '*' by nothing, '2' or '_long' as necessary]
-scaffolding <yes|no>           : scaffolding of contigs used paired end information (default: on)
-max_branch_length <integer>    : maximum length in base pair of bubble (default: 100)
-max_divergence <floating-point>: maximum divergence rate between two branches in a bubble (default: 0.2)
-max_gap_count <integer>        : maximum number of gaps allowed in the alignment of the two branches of a bubble (default: 3)
-min_pair_count <integer>       : minimum number of paired end connections to justify the scaffolding of two long contigs (default: 5)
-max_coverage <floating point>  : removal of high coverage nodes AFTER tour bus (default: no removal)
-coverage_mask <int>    : minimum coverage required for confident regions of contigs (default: 1)
-long_mult_cutoff <int>         : minimum number of long reads required to merge contigs (default: 2)
-unused_reads <yes|no>          : export unused reads in UnusedReads.fa file (default: no)
-alignments <yes|no>            : export a summary of contig alignment to the reference sequences (default: no)
-exportFiltered <yes|no>        : export the long nodes which were eliminated by the coverage filters (default: no)
-clean <yes|no>                 : remove all the intermediary files which are useless for recalculation (default : no)
-very_clean <yes|no>            : remove all the intermediary files (no recalculation possible) (default: no)
-paired_exp_fraction <float>   : remove all the paired end connections which less than the specified fraction of the expected count (default: 0.1)
-shortMatePaired* <yes|no>      : for mate-pair libraries, indicate that the library might be contaminated with paired-end reads (default no)
-conserveLong <yes|no>          : preserve sequences with long reads in them (default no)

Output:
wk_folder/contigs.fa            : fasta file of contigs longer than twice hash length
wk_folder/stats.txt             : stats file (tab-spaced) useful for determining appropriate coverage cutoff
wk_folder/LastGraph             : special formatted file with all the information on the final graph
wk_folder/velvet_asm.afg        : (if requested) AMOS compatible assembly file
    
Example: 
./velvetg wk_folder [options]
string wk_folder; #folder name for files to work on and to save results
float cov_cutoff; #removal of low coverage nodes AFTER tour bus or allow the system to infer it (default: no removal)
int ins_length; #expected distance between two paired end reads (default: no read pairing)
int read_trkg; # (1=yes|0=no) tracking of short read positions in assembly (default:0)
int min_contig_length; #minimum contig length exported to contigs.fa file (default: hash length * 2)
int amos_file; # (1=yes|0=no) #export assembly to AMOS file (default: 0)
float exp_cov; # <floating point|auto>, expected coverage of unique regions or allow the system to infer it (default: no long or paired-end read resolution)
float long_cov_cutoff; #removal of nodes with low long-read coverage AFTER tour bus(default: no removal)


=item Definition

=begin html

<pre>
a reference to a hash where the following keys are defined:
wk_folder has a value which is a string
cov_cutoff has a value which is a float
ins_length has a value which is an int
read_trkg has a value which is an int
min_contig_length has a value which is an int
amos_file has a value which is an int
exp_cov has a value which is a float
long_cov_cutoff has a value which is a float

</pre>

=end html

=begin text

a reference to a hash where the following keys are defined:
wk_folder has a value which is a string
cov_cutoff has a value which is a float
ins_length has a value which is an int
read_trkg has a value which is an int
min_contig_length has a value which is an int
amos_file has a value which is an int
exp_cov has a value which is a float
long_cov_cutoff has a value which is a float


=end text

=back



=head2 VelvetgResults

=over 4



=item Definition

=begin html

<pre>
a reference to a hash where the following keys are defined:
wk_folder has a value which is a string

</pre>

=end html

=begin text

a reference to a hash where the following keys are defined:
wk_folder has a value which is a string


=end text

=back



=cut

package Velvet::VelvetClient::RpcClient;
use base 'JSON::RPC::Client';
use POSIX;
use strict;

#
# Override JSON::RPC::Client::call because it doesn't handle error returns properly.
#

sub call {
    my ($self, $uri, $headers, $obj) = @_;
    my $result;


    {
	if ($uri =~ /\?/) {
	    $result = $self->_get($uri);
	}
	else {
	    Carp::croak "not hashref." unless (ref $obj eq 'HASH');
	    $result = $self->_post($uri, $headers, $obj);
	}

    }

    my $service = $obj->{method} =~ /^system\./ if ( $obj );

    $self->status_line($result->status_line);

    if ($result->is_success) {

        return unless($result->content); # notification?

        if ($service) {
            return JSON::RPC::ServiceObject->new($result, $self->json);
        }

        return JSON::RPC::ReturnObject->new($result, $self->json);
    }
    elsif ($result->content_type eq 'application/json')
    {
        return JSON::RPC::ReturnObject->new($result, $self->json);
    }
    else {
        return;
    }
}


sub _post {
    my ($self, $uri, $headers, $obj) = @_;
    my $json = $self->json;

    $obj->{version} ||= $self->{version} || '1.1';

    if ($obj->{version} eq '1.0') {
        delete $obj->{version};
        if (exists $obj->{id}) {
            $self->id($obj->{id}) if ($obj->{id}); # if undef, it is notification.
        }
        else {
            $obj->{id} = $self->id || ($self->id('JSON::RPC::Client'));
        }
    }
    else {
        # $obj->{id} = $self->id if (defined $self->id);
	# Assign a random number to the id if one hasn't been set
	$obj->{id} = (defined $self->id) ? $self->id : substr(rand(),2);
    }

    my $content = $json->encode($obj);

    $self->ua->post(
        $uri,
        Content_Type   => $self->{content_type},
        Content        => $content,
        Accept         => 'application/json',
	@$headers,
	($self->{token} ? (Authorization => $self->{token}) : ()),
    );
}



1;
