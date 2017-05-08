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


Name of module: Velvet

   This is a KBase module that wraps the open source package "Short read de novo assembler using de Bruijn graphs"
   Velvet_1.2.10

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
$output is an int
VelvethParams is a reference to a hash where the following keys are defined:
	out_folder has a value which is a string
	workspace_name has a value which is a string
	hash_length has a value which is an int
	reads_channels has a value which is a reference to a list where each element is a Velvet.ReadsChannel
ReadsChannel is a reference to a hash where the following keys are defined:
	read_type has a value which is a string
	file_format has a value which is a string
	read_file_info has a value which is a Velvet.ReadFileInfo
	file_layout has a value which is a string
	read_reference has a value which is a Velvet.bool
ReadFileInfo is a reference to a hash where the following keys are defined:
	read_file has a value which is a string
	reference_file has a value which is a string
	left_file has a value which is a string
	right_file has a value which is a string
bool is an int

</pre>

=end html

=begin text

$params is a Velvet.VelvethParams
$output is an int
VelvethParams is a reference to a hash where the following keys are defined:
	out_folder has a value which is a string
	workspace_name has a value which is a string
	hash_length has a value which is an int
	reads_channels has a value which is a reference to a list where each element is a Velvet.ReadsChannel
ReadsChannel is a reference to a hash where the following keys are defined:
	read_type has a value which is a string
	file_format has a value which is a string
	read_file_info has a value which is a Velvet.ReadFileInfo
	file_layout has a value which is a string
	read_reference has a value which is a Velvet.bool
ReadFileInfo is a reference to a hash where the following keys are defined:
	read_file has a value which is a string
	reference_file has a value which is a string
	left_file has a value which is a string
	right_file has a value which is a string
bool is an int


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
$output is an int
VelvetgParams is a reference to a hash where the following keys are defined:
	workspace_name has a value which is a string
	wk_folder has a value which is a string
	output_contigset_name has a value which is a string
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

$params is a Velvet.VelvetgParams
$output is an int
VelvetgParams is a reference to a hash where the following keys are defined:
	workspace_name has a value which is a string
	wk_folder has a value which is a string
	output_contigset_name has a value which is a string
	cov_cutoff has a value which is a float
	ins_length has a value which is an int
	read_trkg has a value which is an int
	min_contig_length has a value which is an int
	amos_file has a value which is an int
	exp_cov has a value which is a float
	long_cov_cutoff has a value which is a float


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
 


=head2 run_velvet

  $output = $obj->run_velvet($params)

=over 4

=item Parameter and return types

=begin html

<pre>
$params is a Velvet.VelvetParams
$output is a Velvet.VelvetResults
VelvetParams is a reference to a hash where the following keys are defined:
	h_params has a value which is a Velvet.VelvethParams
	g_params has a value which is a Velvet.VelvetgParams
VelvethParams is a reference to a hash where the following keys are defined:
	out_folder has a value which is a string
	workspace_name has a value which is a string
	hash_length has a value which is an int
	reads_channels has a value which is a reference to a list where each element is a Velvet.ReadsChannel
ReadsChannel is a reference to a hash where the following keys are defined:
	read_type has a value which is a string
	file_format has a value which is a string
	read_file_info has a value which is a Velvet.ReadFileInfo
	file_layout has a value which is a string
	read_reference has a value which is a Velvet.bool
ReadFileInfo is a reference to a hash where the following keys are defined:
	read_file has a value which is a string
	reference_file has a value which is a string
	left_file has a value which is a string
	right_file has a value which is a string
bool is an int
VelvetgParams is a reference to a hash where the following keys are defined:
	workspace_name has a value which is a string
	wk_folder has a value which is a string
	output_contigset_name has a value which is a string
	cov_cutoff has a value which is a float
	ins_length has a value which is an int
	read_trkg has a value which is an int
	min_contig_length has a value which is an int
	amos_file has a value which is an int
	exp_cov has a value which is a float
	long_cov_cutoff has a value which is a float
VelvetResults is a reference to a hash where the following keys are defined:
	report_name has a value which is a string
	report_ref has a value which is a string

</pre>

=end html

=begin text

$params is a Velvet.VelvetParams
$output is a Velvet.VelvetResults
VelvetParams is a reference to a hash where the following keys are defined:
	h_params has a value which is a Velvet.VelvethParams
	g_params has a value which is a Velvet.VelvetgParams
VelvethParams is a reference to a hash where the following keys are defined:
	out_folder has a value which is a string
	workspace_name has a value which is a string
	hash_length has a value which is an int
	reads_channels has a value which is a reference to a list where each element is a Velvet.ReadsChannel
ReadsChannel is a reference to a hash where the following keys are defined:
	read_type has a value which is a string
	file_format has a value which is a string
	read_file_info has a value which is a Velvet.ReadFileInfo
	file_layout has a value which is a string
	read_reference has a value which is a Velvet.bool
ReadFileInfo is a reference to a hash where the following keys are defined:
	read_file has a value which is a string
	reference_file has a value which is a string
	left_file has a value which is a string
	right_file has a value which is a string
bool is an int
VelvetgParams is a reference to a hash where the following keys are defined:
	workspace_name has a value which is a string
	wk_folder has a value which is a string
	output_contigset_name has a value which is a string
	cov_cutoff has a value which is a float
	ins_length has a value which is an int
	read_trkg has a value which is an int
	min_contig_length has a value which is an int
	amos_file has a value which is an int
	exp_cov has a value which is a float
	long_cov_cutoff has a value which is a float
VelvetResults is a reference to a hash where the following keys are defined:
	report_name has a value which is a string
	report_ref has a value which is a string


=end text

=item Description

Definition of run_velvet

=back

=cut

 sub run_velvet
{
    my($self, @args) = @_;

# Authentication: required

    if ((my $n = @args) != 1)
    {
	Bio::KBase::Exceptions::ArgumentValidationError->throw(error =>
							       "Invalid argument count for function run_velvet (received $n, expecting 1)");
    }
    {
	my($params) = @args;

	my @_bad_arguments;
        (ref($params) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument 1 \"params\" (value was \"$params\")");
        if (@_bad_arguments) {
	    my $msg = "Invalid arguments passed to run_velvet:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	    Bio::KBase::Exceptions::ArgumentValidationError->throw(error => $msg,
								   method_name => 'run_velvet');
	}
    }

    my $url = $self->{url};
    my $result = $self->{client}->call($url, $self->{headers}, {
	    method => "Velvet.run_velvet",
	    params => \@args,
    });
    if ($result) {
	if ($result->is_error) {
	    Bio::KBase::Exceptions::JSONRPC->throw(error => $result->error_message,
					       code => $result->content->{error}->{code},
					       method_name => 'run_velvet',
					       data => $result->content->{error}->{error} # JSON::RPC::ReturnObject only supports JSONRPC 1.1 or 1.O
					      );
	} else {
	    return wantarray ? @{$result->result} : $result->result->[0];
	}
    } else {
        Bio::KBase::Exceptions::HTTP->throw(error => "Error invoking method run_velvet",
					    status_line => $self->{client}->status_line,
					    method_name => 'run_velvet',
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
                method_name => 'run_velvet',
            );
        } else {
            return wantarray ? @{$result->result} : $result->result->[0];
        }
    } else {
        Bio::KBase::Exceptions::HTTP->throw(
            error => "Error invoking method run_velvet",
            status_line => $self->{client}->status_line,
            method_name => 'run_velvet',
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



=head2 bool

=over 4



=item Description

A boolean - 0 for false, 1 for true.
@range (0, 1)


=item Definition

=begin html

<pre>
an int
</pre>

=end html

=begin text

an int

=end text

=back



=head2 ReadFileInfo

=over 4



=item Description

Define a structure that holds the read file name and its use.
Note: only read_file_name is required, the rest are optional.
e.g., 
{"reference_file" => "test_reference.fa", "read_file_name" => "mySortedReads.sam", 
"left_file" => "left.fa", "right_file" => "right.fa"}


=item Definition

=begin html

<pre>
a reference to a hash where the following keys are defined:
read_file has a value which is a string
reference_file has a value which is a string
left_file has a value which is a string
right_file has a value which is a string

</pre>

=end html

=begin text

a reference to a hash where the following keys are defined:
read_file has a value which is a string
reference_file has a value which is a string
left_file has a value which is a string
right_file has a value which is a string


=end text

=back



=head2 ReadsChannel

=over 4



=item Description

Define a structure that mimics the concept of "channel" used by the Velvet program.
string read_type - the read type, e.g., -short, -shortPaired, short2, shortPaired2, -long, or -longPaired
string file_format - the format of the input file, e.g., -fasta, -fastq, -raw,-fasta.gz, -fastq.gz, -raw.gz, -sam, -bam, -fmtAuto
string read_file_info - the hash that holds the details about the read file
string file_layout - the layout of the file, e.g., -interleaved or -separate 
bool read_reference - indicating if a reference file is used


=item Definition

=begin html

<pre>
a reference to a hash where the following keys are defined:
read_type has a value which is a string
file_format has a value which is a string
read_file_info has a value which is a Velvet.ReadFileInfo
file_layout has a value which is a string
read_reference has a value which is a Velvet.bool

</pre>

=end html

=begin text

a reference to a hash where the following keys are defined:
read_type has a value which is a string
file_format has a value which is a string
read_file_info has a value which is a Velvet.ReadFileInfo
file_layout has a value which is a string
read_reference has a value which is a Velvet.bool


=end text

=back



=head2 VelvethParams

=over 4



=item Description

Arguments for velveth input
string workspace_name - the name of the workspace for input/output
string out_folder - the folder name for output files
int hash_length - EITHER an odd integer (if even, it will be decremented) <= 31 (if above, will be reduced)L


=item Definition

=begin html

<pre>
a reference to a hash where the following keys are defined:
out_folder has a value which is a string
workspace_name has a value which is a string
hash_length has a value which is an int
reads_channels has a value which is a reference to a list where each element is a Velvet.ReadsChannel

</pre>

=end html

=begin text

a reference to a hash where the following keys are defined:
out_folder has a value which is a string
workspace_name has a value which is a string
hash_length has a value which is an int
reads_channels has a value which is a reference to a list where each element is a Velvet.ReadsChannel


=end text

=back



=head2 VelvetgParams

=over 4



=item Description

Arguments for run_velvetg

string workspace_name - the name of the workspace from which to take input and store output.
string wk_folder - the name of the folder where the velvet results are created and saved
output_contigset_name - the name of the output contigset list<paired_end_lib>
float cov_cutoff - the removal of low coverage nodes AFTER tour bus or allow the system to infer it (default: no removal)
int ins_length - expected distance between two paired end reads (default: no read pairing)
int read_trkg; -  (1=yes|0=no) tracking of short read positions in assembly (default:0)
int min_contig_length - minimum contig length exported to contigs.fa file (default: hash length * 2)
int amos_file - (1=yes|0=no) #export assembly to AMOS file (default: 0)
float exp_cov - <floating point|auto>, expected coverage of unique regions or allow the system to infer it (default: no long or paired-end read resolution)
float long_cov_cutoff - removal of nodes with low long-read coverage AFTER tour bus(default: no removal)


=item Definition

=begin html

<pre>
a reference to a hash where the following keys are defined:
workspace_name has a value which is a string
wk_folder has a value which is a string
output_contigset_name has a value which is a string
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
workspace_name has a value which is a string
wk_folder has a value which is a string
output_contigset_name has a value which is a string
cov_cutoff has a value which is a float
ins_length has a value which is an int
read_trkg has a value which is an int
min_contig_length has a value which is an int
amos_file has a value which is an int
exp_cov has a value which is a float
long_cov_cutoff has a value which is a float


=end text

=back



=head2 VelvetParams

=over 4



=item Description

Arguments for run_velvet


=item Definition

=begin html

<pre>
a reference to a hash where the following keys are defined:
h_params has a value which is a Velvet.VelvethParams
g_params has a value which is a Velvet.VelvetgParams

</pre>

=end html

=begin text

a reference to a hash where the following keys are defined:
h_params has a value which is a Velvet.VelvethParams
g_params has a value which is a Velvet.VelvetgParams


=end text

=back



=head2 VelvetResults

=over 4



=item Description

Output parameter items for run_velvet

    report_name - the name of the KBaseReport.Report workspace object.
    report_ref - the workspace reference of the report.


=item Definition

=begin html

<pre>
a reference to a hash where the following keys are defined:
report_name has a value which is a string
report_ref has a value which is a string

</pre>

=end html

=begin text

a reference to a hash where the following keys are defined:
report_name has a value which is a string
report_ref has a value which is a string


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
